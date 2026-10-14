"""DuIvy 公众号 Logo 矢量化重制：手绘 JPG -> 亚像素解混 -> 专业 SVG。

用法： python vectorize_logo.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from PIL import Image
from scipy import ndimage
from scipy.spatial import cKDTree

# ------------------------------------------------------------------ 可调参数
SRC = Path(r"C:\Users\hhhhh\Desktop\DuIvyInteraction\logo_duivy\微信图片_20261012002604.jpg")
OUT_DIR = Path(r"C:\Users\hhhhh\Desktop\DuIvyInteraction\logo_duivy")
SVG_OUT = OUT_DIR / "duivy_logo.svg"

CANVAS = 1080          # 输出画布边长，与原始头像规格一致
DPI = 100              # 渲染分辨率，CANVAS / DPI = 英寸数
LEVEL = 0.5            # 覆盖度等值线阈值
MIN_AREA = 800.0       # 小于此面积的连通块视为碎屑丢弃；800 剔除左上红色小尖(323px)，最小真实块为 2645px
SAMPLE_STEP = 0.5      # 轮廓重采样步长（px）
SMOOTH_CORNER = 0.8    # 找角点用的轻平滑 sigma（px）
SMOOTH_FLANK = 3.0     # 角点之间侧翼的重平滑 sigma（px）
FLANK_CAP_DIV = 5.0    # 侧翼 sigma 上限 = 该段弧长 / 此值（保护短段）
FIT_STEP = 1.0         # 拟合前均匀重采样步长（px）
FIT_TOL = 1.40         # 贝塞尔拟合容差（px），取原稿噪声幅度量级
STRAIGHT_TOL = 1.50    # 弦偏小于此值的段吸附为真直线（px）
STRAIGHT_REL = 0.02    # 直线判据的相对项：允许弧长的 2% 作为弦偏
CIRCLE_TOL = 3.50      # 圆拟合最大残差小于此值的段吸附为真圆弧（px）
CIRCLE_MIN_LEN = 80.0  # 启用圆吸附的最短段弧长（px）
ARC_MAX_DEG = 90.0     # 单段圆弧贝塞尔覆盖的最大圆心角（度）
TANGENT_SPAN = 8       # 端点/切分点切向估计基线（采样点，FIT_STEP=1px）
CORNER_CURV = 2.0      # 角点最小局部曲率（度/px），约 R < 28px 才算角点
CORNER_TURN_DEG = 35.0 # 角点在 ±8px 窗口内的最小总转角（度）
CORNER_SPAN_SMALL = 4  # 局部曲率窗口半径（采样点，4 x 0.5px = 2px）
CORNER_SPAN_BIG = 16   # 总转角窗口半径（采样点，16 x 0.5px = 8px）
CORNER_GAP = 32        # 角点合并距离（采样点）
LINE_MERGE_DEG = 30.0  # 直线延伸/合并允许的弱折角（度）
LINE_EXTEND_LEN = 80.0 # 允许被直线延伸吞并的短段弧长上限（px）
SAIL_TIP_OLD = (727.0, 838.0)  # 帆（右侧月牙）下尖位置，用于定位帆轮廓
STERN_POINT = (753.0, 711.0)   # 船尾尖点，用于定位船体轮廓
DECK_TILT_DEG = 4.0            # 船上表面直线绕最左交点的旋转角度（度）
DECK_TILT_FALLOFF = 170.0      # 旋转影响随离该直线距离衰减的尺度（px）
RED_TOUCH_OVERLAP = 2.0        # 红色 T 下端落到船沿时多压入的量（px），保证无缝
RED_DROP_SPAN = 120.0          # 红色 T 下端下移的渐变长度（px）
EXPORT_SIZES = (1080, 512, 256, 128)  # 正式导出的 PNG 边长（px，白底）
SEAM_INSIDE = 0.6              # 提取接缝时要求落在图案内部的墨迹量阈值
SEAM_MATCH = 2.5               # 某段平均距离小于此值即认为它贴着接缝（px）
SEAM_FIT_TOL = 1.4             # 接缝曲线自身的拟合容差（px）；两侧共用同一数据，故不影响缝宽
SEAM_MIN_LEN = 20.0            # 接缝折线的最短有效长度（px）
HANDLE_MAX = 3.0       # 手柄长度上限（相对弦长），防止病态解跑飞
NEWTON_ITERS = 4       # 最小二乘重参数化迭代次数
MAX_DEPTH = 12         # 递归细分深度上限
DECIMALS = 2           # SVG 坐标保留小数位
CLASS_TOL = 40.0       # 渲染图按色距分类的容差
KINK_DEG = 5.0         # 判为折角的切向夹角阈值（度）

WHITE = (255.0, 255.0, 255.0)
BLUE_HEX = "#1F6FD2"   # 自原稿内部像素实测
RED_HEX = "#FC6671"    # 自原稿内部像素实测

SVG_HEAD = (
    '<svg xmlns="http://www.w3.org/2000/svg" width="{c}" height="{c}" '
    'viewBox="0 0 {c} {c}" shape-rendering="geometricPrecision">\n'
    "<title>DuIvy 杜艾维</title>\n"
    "<desc>字母 D 与粉色 T 构成的帆船意象 Logo</desc>\n"
)


# ------------------------------------------------------------------ 基础工具
def rgb_to_float(hex_color: str) -> np.ndarray:
    """'#RRGGBB' -> (3,) float RGB。"""
    h = hex_color.lstrip("#")
    return np.array([int(h[i : i + 2], 16) for i in (0, 2, 4)], dtype=float)


def load_rgb(path: Path) -> np.ndarray:
    """读取图片为 (H, W, 3) float RGB。"""
    return np.asarray(Image.open(path).convert("RGB"), dtype=float)


def unmix(rgb: np.ndarray, blue: np.ndarray, red: np.ndarray, white: np.ndarray):
    """把抗锯齿图像解混为蓝、红两张亚像素覆盖度场。"""
    basis = np.stack([blue - white, red - white])
    gram = basis @ basis.T
    rhs = (rgb - white) @ basis.T
    coef = np.linalg.solve(gram, rhs[..., None])[..., 0]
    return np.clip(coef[..., 0], 0.0, 1.0), np.clip(coef[..., 1], 0.0, 1.0)


def seam_fields(alpha_blue: np.ndarray, alpha_red: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """重构红蓝场，使两者在交界处共享同一条接缝（无白缝、无重叠）。"""
    weight = np.clip(alpha_blue + alpha_red, 0.0, 1.0)
    diff = alpha_blue - alpha_red
    blue = np.clip(np.minimum(weight, 0.5 + diff), 0.0, 1.0)
    red = np.clip(np.minimum(weight, 0.5 - diff), 0.0, 1.0)
    return blue, red


def drop_specks(alpha: np.ndarray, level: float, min_area: float) -> np.ndarray:
    """抹掉面积小于 min_area 的连通块。"""
    binary = alpha > level
    lab, n = ndimage.label(binary)
    if n == 0:
        return alpha
    sizes = ndimage.sum(binary, lab, range(1, n + 1))
    keep = np.isin(lab, np.nonzero(sizes >= min_area)[0] + 1)
    return np.where(keep, alpha, 0.0)


def contour_loops(alpha: np.ndarray, level: float) -> list[np.ndarray]:
    """提取等值线的亚像素闭合轮廓。"""
    fig = Figure()
    FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    segs = ax.contour(np.pad(alpha, 1), levels=[level]).allsegs[0]
    loops = []
    for seg in segs:
        pts = np.asarray(seg, dtype=float) - 1.0
        if len(pts) > 1 and np.allclose(pts[0], pts[-1]):
            pts = pts[:-1]
        if len(pts) >= 4:
            loops.append(pts)
    return loops


# ------------------------------------------------------------------ 轮廓处理
def resample_closed(pts: np.ndarray, step: float) -> np.ndarray:
    """按固定弧长重采样闭合轮廓。"""
    ring = np.vstack([pts, pts[:1]])
    seg = np.hypot(*np.diff(ring, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 0.0:
        return pts
    t = np.linspace(0.0, s[-1], max(8, int(round(s[-1] / step))), endpoint=False)
    return np.stack([np.interp(t, s, ring[:, 0]), np.interp(t, s, ring[:, 1])], axis=1)


def resample_open(pts: np.ndarray, step: float) -> np.ndarray:
    """按固定弧长重采样开放折线，使拟合点距均匀。"""
    seg = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 0.0:
        return pts
    t = np.linspace(0.0, s[-1], max(3, int(round(s[-1] / step)) + 1))
    return np.stack([np.interp(t, s, pts[:, 0]), np.interp(t, s, pts[:, 1])], axis=1)


def smooth_closed(pts: np.ndarray, step: float, sigma_px: float) -> np.ndarray:
    """沿弧长做循环高斯平滑，去掉锯齿抖动。"""
    ring = resample_closed(pts, step)
    sigma = max(sigma_px / step, 1e-6)
    return ndimage.gaussian_filter1d(ring, sigma=sigma, axis=0, mode="wrap")


def turn_cosines(pts: np.ndarray, span: int) -> np.ndarray:
    """每个顶点在正负 span 点窗口内的转角余弦。"""
    prev = pts - np.roll(pts, span, axis=0)
    nxt = np.roll(pts, -span, axis=0) - pts
    denom = np.hypot(*prev.T) * np.hypot(*nxt.T)
    cos = np.ones(len(pts))
    ok = denom > 0.0
    cos[ok] = np.einsum("ij,ij->i", prev[ok], nxt[ok]) / denom[ok]
    return cos


def corners_of(pts: np.ndarray, gap: int) -> list[int]:
    """按局部曲率与总转角判定角点，并合并彼此过近的候选。"""
    small = np.degrees(np.arccos(np.clip(turn_cosines(pts, CORNER_SPAN_SMALL), -1.0, 1.0)))
    big = np.degrees(np.arccos(np.clip(turn_cosines(pts, CORNER_SPAN_BIG), -1.0, 1.0)))
    curv = small / (2.0 * CORNER_SPAN_SMALL * SAMPLE_STEP)
    cand = np.nonzero((curv > CORNER_CURV) & (big > CORNER_TURN_DEG))[0]
    kept: list[int] = []
    for i in sorted((int(j) for j in cand), key=lambda j: -curv[j]):
        if all(min((i - j) % len(curv), (j - i) % len(curv)) >= gap for j in kept):
            kept.append(i)
    return sorted(kept)


def split_runs(pts: np.ndarray, corners: list[int]) -> list[np.ndarray]:
    """按角点把闭合轮廓切成若干开放段。"""
    idx = sorted(set(corners))
    runs = []
    for j, start in enumerate(idx):
        end = idx[(j + 1) % len(idx)]
        if end > start:
            run = pts[start : end + 1]
        else:
            run = np.vstack([pts[start:], pts[: end + 1]])
        if len(run) >= 2:
            runs.append(run)
    return runs


# ------------------------------------------------------------------ 曲线拟合
def _chord_u(pts: np.ndarray) -> np.ndarray:
    """弦长参数化。"""
    d = np.hypot(*np.diff(pts, axis=0).T)
    u = np.concatenate([[0.0], np.cumsum(d)])
    return u / u[-1] if u[-1] > 0.0 else np.linspace(0.0, 1.0, len(pts))


def bezier_eval(bez: np.ndarray, u: np.ndarray) -> np.ndarray:
    """求三次贝塞尔曲线在参数 u 处的点。"""
    c = 3.0 * (bez[1] - bez[0])
    d = 3.0 * (bez[2] - bez[1]) - c
    e = bez[3] - bez[0] - c - d
    t = u[:, None]
    return bez[0] + c * t + d * t * t + e * t * t * t


def fit_segment_once(pts: np.ndarray, u: np.ndarray, t1: np.ndarray, t2: np.ndarray) -> np.ndarray:
    """最小二乘求单段三次贝塞尔的控制点（Schneider）。"""
    b0, b1 = (1 - u) ** 3, 3 * (1 - u) ** 2 * u
    b2, b3 = 3 * (1 - u) * u**2, u**3
    p0, p3 = pts[0], pts[-1]
    chord = float(np.hypot(*(p3 - p0))) or 1e-9
    tmp = pts - (p0 * (b0 + b1)[:, None] + p3 * (b2 + b3)[:, None])
    c00 = float(t1 @ t1) * float(b1 @ b1)
    c01 = float(t1 @ t2) * float(b1 @ b2)
    c11 = float(t2 @ t2) * float(b2 @ b2)
    x0, x1 = float(b1 @ (tmp @ t1)), float(b2 @ (tmp @ t2))
    det = c00 * c11 - c01 * c01
    a_l = a_r = -1.0
    if abs(det) > 1e-12:
        a_l = (x0 * c11 - x1 * c01) / det
        a_r = (x1 * c00 - x0 * c01) / det
    limit = HANDLE_MAX * chord
    if not 0.0 < a_l <= limit or not 0.0 < a_r <= limit:
        a_l = a_r = chord / 3.0
    return np.stack([p0, p0 + t1 * a_l, p3 + t2 * a_r, p3])


def max_error(bez: np.ndarray, pts: np.ndarray, u: np.ndarray) -> tuple[float, int]:
    """返回最大偏差及其位置下标。"""
    d = np.hypot(*(bezier_eval(bez, u) - pts).T)
    i = int(np.argmax(d))
    return float(d[i]), i


def reparameterize(bez: np.ndarray, pts: np.ndarray, u: np.ndarray) -> np.ndarray:
    """牛顿迭代重参数化。"""
    for _ in range(NEWTON_ITERS):
        q = bezier_eval(bez, u)
        c = 3.0 * (bez[1] - bez[0])
        d = 3.0 * (bez[2] - bez[1]) - c
        e = bez[3] - bez[0] - c - d
        t = u[:, None]
        q1 = c + 2 * d * t + 3 * e * t * t
        q2 = 2 * d + 6 * e * t
        num = np.einsum("ij,ij->i", q - pts, q1)
        den = np.einsum("ij,ij->i", q1, q1) + np.einsum("ij,ij->i", q - pts, q2)
        u = np.clip(u - num / np.where(np.abs(den) < 1e-12, 1e-12, den), 0.0, 1.0)
    return u


def fit_curve(pts: np.ndarray, t1: np.ndarray, t2: np.ndarray, tol: float, depth: int = 0):
    """递归拟合一段点列为若干三次贝塞尔段。"""
    if len(pts) == 2:
        dist = float(np.hypot(*(pts[1] - pts[0]))) / 3.0
        return [np.stack([pts[0], pts[0] + t1 * dist, pts[1] + t2 * dist, pts[1]])]
    u = _chord_u(pts)
    bez = fit_segment_once(pts, u, t1, t2)
    err, i = max_error(bez, pts, u)
    if err < tol:
        return [bez]
    if err < tol * tol or depth >= MAX_DEPTH:
        for _ in range(NEWTON_ITERS):
            u = reparameterize(bez, pts, u)
            bez = fit_segment_once(pts, u, t1, t2)
        err, i = max_error(bez, pts, u)
        if err < tol or depth >= MAX_DEPTH:
            return [bez]
    i = min(max(i, 1), len(pts) - 2)
    k = max(1, min(TANGENT_SPAN, i, len(pts) - 1 - i))
    mid = pts[i - k] - pts[i + k]
    nrm = float(np.hypot(*mid))
    mid = mid / nrm if nrm > 0.0 else np.array([1.0, 0.0])
    return fit_curve(pts[: i + 1], t1, mid, tol, depth + 1) + fit_curve(
        pts[i:], -mid, t2, tol, depth + 1
    )


def _unit(v: np.ndarray, fallback: np.ndarray) -> np.ndarray:
    """单位化，零向量时退回 fallback。"""
    nrm = float(np.hypot(*v))
    return v / nrm if nrm > 1e-12 else fallback


def smooth_open(pts: np.ndarray, step: float, sigma_px: float) -> np.ndarray:
    """按固定弧长重采样并平滑开放折线，两端点保持不变。"""
    rs = resample_open(pts, step)
    out = ndimage.gaussian_filter1d(rs, sigma=max(sigma_px / step, 1e-6), axis=0, mode="nearest")
    out[0], out[-1] = pts[0], pts[-1]
    return out


def _arc_length(pts: np.ndarray) -> float:
    """折线总长。"""
    return float(np.hypot(*np.diff(pts, axis=0).T).sum())


def straight_segment(p0: np.ndarray, p1: np.ndarray) -> np.ndarray:
    """两点间的直线，以三次贝塞尔表示。"""
    return np.stack([p0, p0 + (p1 - p0) / 3.0, p0 + 2.0 * (p1 - p0) / 3.0, p1])


def max_chord_dev(pts: np.ndarray) -> float:
    """各点到首尾弦的最大距离。"""
    a, ab = pts[0], pts[-1] - pts[0]
    nrm = float(np.hypot(*ab))
    if nrm < 1e-9:
        return float(np.hypot(*(pts - a).T).max())
    return float(np.abs(ab[0] * (pts[:, 1] - a[1]) - ab[1] * (pts[:, 0] - a[0])).max() / nrm)


def fit_circle(pts: np.ndarray) -> tuple[np.ndarray, float, float]:
    """最小二乘拟合圆，返回圆心、半径与最大残差。"""
    mat = np.column_stack([2.0 * pts[:, 0], 2.0 * pts[:, 1], np.ones(len(pts))])
    sol, *_ = np.linalg.lstsq(mat, (pts**2).sum(axis=1), rcond=None)
    center = sol[:2]
    radius = float(np.sqrt(max(sol[2] + float(center @ center), 0.0)))
    res = float(np.abs(np.hypot(*(pts - center).T) - radius).max())
    return center, radius, res


def force_ends(segs: list[np.ndarray], p0: np.ndarray, p1: np.ndarray) -> list[np.ndarray]:
    """把首末端点精确对齐到给定点，只平移端部手柄。"""
    out = [seg.copy() for seg in segs]
    out[0][1] += p0 - out[0][0]
    out[0][0] = p0
    out[-1][2] += p1 - out[-1][3]
    out[-1][3] = p1
    return out


def arc_beziers(center: np.ndarray, radius: float, a0: float, a1: float) -> list[np.ndarray]:
    """按圆心角把圆弧拆成若干三次贝塞尔段。"""
    n = max(1, int(np.ceil(abs(np.degrees(a1 - a0)) / ARC_MAX_DEG)))
    out = []
    for k in range(n):
        t0 = a0 + (a1 - a0) * k / n
        t1 = a0 + (a1 - a0) * (k + 1) / n
        handle = 4.0 / 3.0 * np.tan((t1 - t0) / 4.0) * radius
        p0 = center + radius * np.array([np.cos(t0), np.sin(t0)])
        p3 = center + radius * np.array([np.cos(t1), np.sin(t1)])
        d0 = np.array([-np.sin(t0), np.cos(t0)])
        d1 = np.array([-np.sin(t1), np.cos(t1)])
        out.append(np.stack([p0, p0 + handle * d0, p3 - handle * d1, p3]))
    return out


def circle_arc(pts: np.ndarray) -> list[np.ndarray] | None:
    """该段近似圆弧时返回真圆弧贝塞尔，否则 None。"""
    if _arc_length(pts) < CIRCLE_MIN_LEN:
        return None
    center, radius, res = fit_circle(pts)
    if res >= CIRCLE_TOL or radius < 1.0:
        return None
    ang = np.unwrap(np.arctan2(pts[:, 1] - center[1], pts[:, 0] - center[0]))
    if abs(ang[-1] - ang[0]) > 2.0 * np.pi - 1e-3:
        return None
    return force_ends(arc_beziers(center, radius, ang[0], ang[-1]), pts[0], pts[-1])


def fit_run(pts: np.ndarray, tol: float, snap: bool = True) -> list[np.ndarray]:
    """把一条开放折线拟合成贝塞尔段，近直线吸附直线、近圆弧吸附圆弧。"""
    if len(pts) < 3:
        return []
    if float(np.hypot(*(pts[-1] - pts[0]))) < SAMPLE_STEP:
        mid = len(pts) // 2
        return fit_run(pts[: mid + 1], tol, snap) + fit_run(pts[mid:], tol, snap)
    if snap and max_chord_dev(pts) < max(STRAIGHT_TOL, _arc_length(pts) * STRAIGHT_REL):
        return [straight_segment(pts[0], pts[-1])]
    if snap:
        arc = circle_arc(pts)
        if arc is not None:
            return arc
    k = min(TANGENT_SPAN, len(pts) - 1)
    t1 = _unit(pts[k] - pts[0], np.array([1.0, 0.0]))
    t2 = _unit(pts[-1 - k] - pts[-1], np.array([-1.0, 0.0]))
    return fit_curve(pts, t1, t2, tol)


def _is_line(run: list[np.ndarray]) -> bool:
    """判断该段是否为一条真直线段。"""
    return len(run) == 1 and abs(float(np.cross(run[0][1] - run[0][0], run[0][3] - run[0][0]))) < 1e-6


def _join_deg(a: list[np.ndarray], b: list[np.ndarray]) -> float:
    """相邻两段接缝处的切向夹角（度）。"""
    return _join_angle(a[-1][3] - a[-1][2], b[0][1] - b[0][0])


def _short(run: list[np.ndarray]) -> bool:
    """该段是否短于可被直线延伸吞并的长度。"""
    return _arc_length(run_polyline(run)) < LINE_EXTEND_LEN


def merge_weak_lines(loop: list) -> list:
    """把真直线段跨过弱折角延伸或合并，消除台阶式小折。"""
    runs = list(loop)
    changed = True
    while changed and len(runs) > 2:
        changed = False
        for i in range(len(runs)):
            a, b = runs[i - 1], runs[i]
            if _join_deg(a, b) >= LINE_MERGE_DEG:
                continue
            if _is_line(a) and (_is_line(b) or _short(b)) or _is_line(b) and _short(a):
                runs[i - 1] = [straight_segment(a[0][0], b[-1][-1])]
                runs.pop(i)
                changed = True
                break
    return runs


def _has_end_near(loop: list, target: np.ndarray, tol: float = 12.0) -> bool:
    """该闭合轮廓是否有端点落在给定点附近。"""
    d0 = [float(np.hypot(*(run[0][0] - target))) for run in loop]
    d1 = [float(np.hypot(*(run[-1][-1] - target))) for run in loop]
    return min(d0 + d1) < tol


def _line_dist(pts: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """各点到直线 ab 的垂直距离。"""
    ab = b - a
    normal = np.array([-ab[1], ab[0]]) / np.hypot(*ab)
    return np.abs((pts - a) @ normal)


def _min_line_dist(pts: np.ndarray, segs: list) -> np.ndarray:
    """到若干直线的最小距离。"""
    return np.minimum.reduce([_line_dist(pts, a, b) for a, b in segs])


def _rotate_about(pts: np.ndarray, pivot: np.ndarray, ang: float) -> np.ndarray:
    """把点列绕 pivot 旋转 ang 弧度。"""
    c, s = np.cos(ang), np.sin(ang)
    v = pts - pivot
    return pivot + np.stack([v[:, 0] * c - v[:, 1] * s, v[:, 0] * s + v[:, 1] * c], axis=1)


def _second_intersection(c1: np.ndarray, c2: np.ndarray, known: np.ndarray) -> np.ndarray:
    """两圆已交于 known 点，求另一个交点。"""
    u = np.array([-(c2 - c1)[1], (c2 - c1)[0]])
    nrm = float(np.hypot(*u))
    if nrm < 1e-9:
        return known
    u = u / nrm
    return known - 2.0 * float((known - c1) @ u) * u


def _rearc(run: list, center: np.ndarray, radius: float, p0: np.ndarray, p1: np.ndarray) -> list:
    """把一段曲线重建为指定圆上从 p0 到 p1 的圆弧段。"""
    pl = run_polyline(run)
    ang = np.unwrap(np.arctan2(pl[:, 1] - center[1], pl[:, 0] - center[0]))
    return force_ends(arc_beziers(center, radius, ang[0], ang[-1]), p0, p1)


def _deck_pivot(loops: list, segs: list, sail_tip: np.ndarray) -> np.ndarray:
    """甲板直线与帆内弧所在圆的最左交点，作为旋转支点。"""
    sail = next(loop for loop in loops if _has_end_near(loop, sail_tip))
    inner_pl = min((run_polyline(run) for run in sail), key=lambda pl: float(_min_line_dist(pl, segs).min()))
    center, radius, _ = fit_circle(inner_pl)
    a, b = min(segs, key=lambda s: float(_min_line_dist(inner_pl, [s]).min()))
    d = (b - a) / np.hypot(*(b - a))
    va = a - center
    disc = (2.0 * float(va @ d)) ** 2 - 4.0 * (float(va @ va) - radius * radius)
    if disc <= 0.0:
        return a
    hits = [a + t * d for t in ((-2.0 * float(va @ d) - np.sqrt(disc)) / 2.0, (-2.0 * float(va @ d) + np.sqrt(disc)) / 2.0)]
    return min(hits, key=lambda p: p[0])


def _line_at_x(a: np.ndarray, b: np.ndarray, x: float) -> float:
    """直线 ab 在给定 x 处的 y。"""
    return float(a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0]))


def hull_index(loops: list, stern: np.ndarray) -> int:
    """定位含船尾尖点的轮廓下标。"""
    return next(i for i, loop in enumerate(loops) if _has_end_near(loop, stern))


def deck_segments(hull: list) -> list:
    """取船体轮廓中的甲板直线段（起点，终点）。"""
    return [(run[0][0], run[-1][-1]) for run in hull if _is_line(run)]


def deck_y_at(segs: list, x: float) -> float:
    """甲板线在给定 x 处的高度（取 x 落在其跨度内的那一段）。"""
    seg = min(segs, key=lambda s: max(min(s[0][0], s[1][0]) - x, x - max(s[0][0], s[1][0]), 0.0))
    return _line_at_x(seg[0], seg[1], x)


def drop_red_foot(red_loops: list, bottom_y: float, drop: float, span: float) -> list:
    """把红色轮廓下端整体下移 drop，渐变过渡到上段保持不动。"""
    out = []
    for loop in red_loops:
        runs = []
        for run in loop:
            pts = resample_open(run_polyline(run), FIT_STEP)
            weight = np.clip((pts[:, 1] - (bottom_y - span)) / span, 0.0, 1.0)
            runs.append(fit_run(pts + np.array([0.0, drop]) * weight[:, None], FIT_TOL))
        out.append(runs)
    return out


def red_foot_point(loops: list) -> np.ndarray:
    """红色轮廓上最低的端点（即 T 的下切口处）。"""
    pts = [run[k][0] for loop in loops for run in loop for k in (0, -1)]
    return max(pts, key=lambda p: p[1])


def _nearest_seam(run: list, seams: list) -> tuple[int | None, float]:
    """该段最贴的接缝曲线下标与平均距离。"""
    pl = run_polyline(run)
    best, best_d = None, np.inf
    for i, seam in enumerate(seams):
        d = float(
            np.minimum.reduce([_pt_seg_dist(pl, seam[k], seam[k + 1]) for k in range(len(seam) - 1)]).mean()
        )
        if d < best_d:
            best, best_d = i, d
    return best, best_d


def share_seams(layers: dict, seams: list, match: float, seam_tol: float) -> dict:
    """两色的接缝边统一取同一条接缝曲线的拟合结果：该曲线只拟合一次，两侧共用同一份数据。"""
    fits = [fit_run(resample_open(seam, FIT_STEP), seam_tol) for seam in seams]
    out = {}
    for key, loops in layers.items():
        new_loops = []
        for loop in loops:
            runs = []
            for run in loop:
                idx, mean = _nearest_seam(run, seams)
                if idx is None or mean >= match:
                    runs.append(run)
                    continue
                segs = [seg.copy() for seg in fits[idx]]
                if float(np.hypot(*(run[0][0] - segs[0][0]))) > float(np.hypot(*(run[0][0] - segs[-1][-1]))):
                    segs = reverse_run(segs)
                runs.append(force_ends(segs, run[0][0], run[-1][-1]))
            new_loops.append(runs)
        out[key] = new_loops
    return out


def _cyclic_runs(pts: np.ndarray, mask: np.ndarray) -> list[np.ndarray]:
    """取闭合折线上 mask 连续为真的段落。"""
    if mask.all():
        return [pts]
    start = int(np.argmin(mask))
    runs, cur = [], []
    for k in np.roll(np.arange(len(mask)), -start):
        if mask[k]:
            cur.append(k)
        elif cur:
            runs.append(pts[np.array(cur)])
            cur = []
    if cur:
        runs.append(pts[np.array(cur)])
    return [r for r in runs if len(r) > 1]


def seam_polylines(fields: dict) -> list[np.ndarray]:
    """从解混场提取红蓝分界线（接缝）折线，只保留落在图案内部的段落。"""
    diff = fields["blue"] - fields["red"]
    weight = np.clip(fields["blue"] + fields["red"], 0.0, 1.0)
    out = []
    for poly in contour_loops(diff, 0.0):
        idx = np.clip(poly.astype(int), 0, CANVAS - 1)
        for run in _cyclic_runs(poly, weight[idx[:, 1], idx[:, 0]] > SEAM_INSIDE):
            if _arc_length(run) > SEAM_MIN_LEN:
                out.append(run)
    return out


def tilt_deck(loops: list, hull_idx: int, sail_tip: np.ndarray, angle_deg: float, falloff: float) -> list:
    """把船的上表面直线绕它与帆内弧的最左交点旋转，影响随离该线距离衰减。"""
    hull = loops[hull_idx]
    segs = deck_segments(hull)
    pivot = _deck_pivot(loops, segs, sail_tip)
    ang = np.radians(angle_deg)
    out = []
    for loop in loops:
        if loop is not hull:
            out.append(loop)
            continue
        runs = []
        for run in loop:
            pts = resample_open(run_polyline(run), FIT_STEP)
            u = np.clip(1.0 - _min_line_dist(pts, segs) / falloff, 0.0, 1.0)
            weight = u * u * (3.0 - 2.0 * u)
            runs.append(fit_run(pts + weight[:, None] * (_rotate_about(pts, pivot, ang) - pts), FIT_TOL))
        out.append(runs)
    return out


def _flank_sigma(run: np.ndarray) -> float:
    """按段长限制侧翼平滑强度，避免短段被压平。"""
    return min(SMOOTH_FLANK, _arc_length(run) / FLANK_CAP_DIV)


def vectorize(alpha: np.ndarray) -> list[list[list[np.ndarray]]]:
    """覆盖度场 -> 每条轮廓 -> 每段开放曲线 -> 贝塞尔段。"""
    loops = []
    for contour in contour_loops(drop_specks(alpha, LEVEL, MIN_AREA), LEVEL):
        dense = smooth_closed(contour, SAMPLE_STEP, SMOOTH_CORNER)
        corners = corners_of(dense, CORNER_GAP) or [0]
        runs = []
        for run in split_runs(dense, corners):
            smoothed = smooth_open(run, SAMPLE_STEP, _flank_sigma(run))
            segs = fit_run(resample_open(smoothed, FIT_STEP), FIT_TOL)
            if segs:
                runs.append(segs)
        if runs:
            loops.append(runs)
    return loops


def run_polyline(run: list[np.ndarray], per_seg: int = 16) -> np.ndarray:
    """把一段贝塞尔序列密集采样成折线。"""
    pts = []
    for seg in run:
        u = np.linspace(0.0, 1.0, per_seg, endpoint=False)
        pts.append(bezier_eval(seg, u))
    pts.append(run[-1][-1])
    return np.vstack(pts)


def _pt_seg_dist(pts: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """一组点到单条线段的距离。"""
    ab = b - a
    t = np.clip((pts - a) @ ab / max(float(ab @ ab), 1e-12), 0.0, 1.0)
    proj = a + t[:, None] * ab
    return np.hypot(*(proj - pts).T)


def reverse_run(run: list[np.ndarray]) -> list[np.ndarray]:
    """反转一段贝塞尔序列的方向。"""
    return [np.flip(seg, axis=0) for seg in reversed(run)]


# ------------------------------------------------------------------ 输出与校验
def all_points(loops) -> np.ndarray:
    """收集所有贝塞尔控制点。"""
    return np.vstack([seg for loop in loops for run in loop for seg in run])


def centering_shift(layers: dict, canvas: int) -> np.ndarray:
    """求把图案居中到画布所需的平移量。"""
    pts = np.vstack([all_points(v) for v in layers.values()])
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    return (canvas - (hi - lo)) / 2.0 - lo


def _num(v: float) -> str:
    """格式化单个坐标值。"""
    return f"{v:.{DECIMALS}f}"


def loop_to_d(runs: list[list[np.ndarray]], shift: np.ndarray) -> str:
    """一条闭合轮廓的全部曲线段 -> SVG 的一个闭合子路径。"""
    head = runs[0][0][0] + shift
    parts = [f"M{_num(head[0])} {_num(head[1])}"]
    for run in runs:
        for seg in run:
            p = seg + shift
            parts.append("C" + " ".join(f"{_num(c[0])} {_num(c[1])}" for c in p[1:]))
    parts.append("Z")
    return "".join(parts)


def build_svg(layers: dict, shift: np.ndarray) -> str:
    """拼装最终 SVG。"""
    body = []
    for key, color in (("blue", BLUE_HEX), ("red", RED_HEX)):
        d = "".join(loop_to_d(loop, shift) for loop in layers[key])
        body.append(f'<path fill="{color}" d="{d}"/>')
    return SVG_HEAD.format(c=CANVAS) + "\n".join(body) + "\n</svg>\n"


def loop_to_mpl(runs: list[list[np.ndarray]], shift: np.ndarray) -> MplPath:
    """一条闭合轮廓 -> matplotlib 路径。"""
    start = runs[0][0][0] + shift
    verts, codes = [start], [MplPath.MOVETO]
    for run in runs:
        for seg in run:
            p = seg + shift
            verts.extend([p[1], p[2], p[3]])
            codes.extend([MplPath.CURVE4] * 3)
    verts.append(start)
    codes.append(MplPath.CLOSEPOLY)
    return MplPath(np.asarray(verts), np.asarray(codes))


def _region_axes(background: str, box: tuple[int, int, int, int], scale: int) -> tuple[Figure, Axes]:
    """建一块只覆盖指定区域的离屏绘图区。"""
    x0, y0, x1, y1 = box
    fig = Figure(figsize=((x1 - x0) * scale / DPI, (y1 - y0) * scale / DPI), dpi=DPI)
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.axis("off")
    fig.patch.set_facecolor(background)
    return fig, ax


def _canvas_axes(background: str) -> tuple[Figure, Axes]:
    """建一块与输出画布等像素的离屏绘图区。"""
    return _region_axes(background, (0, 0, CANVAS, CANVAS), 1)


def render_rgba(groups: list[tuple[list, str]], shift: np.ndarray, background: str) -> np.ndarray:
    """按矢量几何渲染成 RGBA（background='none' 即透明底）。"""
    fig, ax = _canvas_axes(background)
    for loops, color in groups:
        for loop in loops:
            ax.add_patch(PathPatch(loop_to_mpl(loop, shift), facecolor=color, edgecolor="none"))
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())


def render_rgb(groups: list[tuple[list, str]], shift: np.ndarray, background: str) -> np.ndarray:
    """把若干 (轮廓组, 颜色) 按矢量几何渲染成抗锯齿图。"""
    return render_rgba(groups, shift, background)[..., :3].astype(float)


def render_coverage(loops, shift: np.ndarray) -> np.ndarray:
    """把矢量几何渲染成抗锯齿覆盖度图。"""
    return render_rgb([(loops, "white")], shift, "black")[..., 0] / 255.0


def write_exports(groups: list[tuple[list, str]], shift: np.ndarray) -> None:
    """导出正式 PNG：白底多尺寸 + 透明底 1080。"""
    final = render_rgb(groups, shift, "white")
    for size in EXPORT_SIZES:
        img = Image.fromarray(np.clip(final, 0, 255).astype(np.uint8))
        if size != CANVAS:
            img = img.resize((size, size), Image.LANCZOS)
        img.save(OUT_DIR / f"duivy_logo_{size}.png")
    Image.fromarray(render_rgba(groups, shift, "none")).save(OUT_DIR / "duivy_logo_transparent_1080.png")


def mask_metrics(ref: np.ndarray, got: np.ndarray) -> tuple[float, float, float, float]:
    """比较两个覆盖度场的形状一致性。"""
    a, b = ref > LEVEL, got > LEVEL
    union = int(np.logical_or(a, b).sum())
    iou = int(np.logical_and(a, b).sum()) / union if union else 1.0
    edge = b & ~ndimage.binary_erosion(b, np.ones((3, 3)))
    dev = ndimage.distance_transform_edt(~a)[edge] if edge.any() else np.zeros(1)
    return iou, float(dev.mean()), float(np.percentile(dev, 95)), float(dev.max())


# ------------------------------------------------------------------ 主流程
def _join_angle(a_out: np.ndarray, a_in: np.ndarray) -> float:
    """两个切向量之间的夹角（度）。"""
    n = float(np.hypot(*a_out)) * float(np.hypot(*a_in))
    if n <= 1e-12:
        return 0.0
    return float(np.degrees(np.arccos(np.clip(float(a_out @ a_in) / n, -1.0, 1.0))))


def within_run_angles(loops: list) -> list[float]:
    """同一段曲线内部各接缝的切向夹角，应接近 0（G1 连续）。"""
    out = []
    for loop in loops:
        for run in loop:
            for i in range(len(run) - 1):
                out.append(_join_angle(run[i][3] - run[i][2], run[i + 1][1] - run[i + 1][0]))
    return out


def loop_kink_angles(loop: list) -> list[float]:
    """一条闭合轮廓所有相邻段接缝的切向夹角。"""
    segs = [seg for run in loop for seg in run]
    return [
        _join_angle(segs[i][3] - segs[i][2], segs[(i + 1) % len(segs)][1] - segs[(i + 1) % len(segs)][0])
        for i in range(len(segs))
    ]


def inflections(loop: list, samples: int = 12) -> int:
    """闭合轮廓上曲率变号次数，用来度量摆动。"""
    pts = np.vstack([run_polyline(run, samples) for run in loop])
    d1 = np.roll(pts, -1, axis=0) - pts
    d2 = np.roll(d1, -1, axis=0) - d1
    sign = np.sign(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0])
    sign = sign[sign != 0.0]
    return int((np.diff(sign) != 0.0).sum())


def smooth_report(name: str, loops: list) -> None:
    """打印段数、段内 G1 夹角、折角数与摆动次数。"""
    segs = sum(len(run) for loop in loops for run in loop)
    runs = sum(len(loop) for loop in loops)
    inner = within_run_angles(loops)
    kinks = [a for loop in loops for a in loop_kink_angles(loop) if a > KINK_DEG]
    wobble = sum(inflections(loop) for loop in loops)
    print(
        f"[{name}] runs={runs} beziers={segs} 段内夹角 max={max(inner) if inner else 0:.2f}° "
        f"折角(>{KINK_DEG:.0f}°)={len(kinks)} 摆动={wobble}"
    )


def kept_blocks(alpha: np.ndarray) -> list[tuple[int, tuple[int, int, int, int]]]:
    """列出面积达标的连通块及其外接框。"""
    binary = alpha > LEVEL
    lab, n = ndimage.label(binary)
    out = []
    for i, size in enumerate(ndimage.sum(binary, lab, range(1, n + 1)), start=1):
        if size >= MIN_AREA:
            yy, xx = np.nonzero(lab == i)
            out.append((int(size), (int(xx.min()), int(yy.min()), int(xx.max()), int(yy.max()))))
    return sorted(out, reverse=True)


def color_stats(rgb: np.ndarray, alpha: np.ndarray) -> tuple[str, str]:
    """内部像素的中位数色与众数色。"""
    inner = rgb[alpha > 0.95]
    med = np.median(inner, axis=0).astype(int)
    q = inner.astype(int) // 4 * 4 + 2
    uniq, cnt = np.unique(q, axis=0, return_counts=True)
    mode = uniq[int(np.argmax(cnt))]
    return "#%02X%02X%02X" % tuple(med), "#%02X%02X%02X" % tuple(mode)


def sliver_px(layers: dict) -> int:
    """统计轮廓内既不属于蓝也不属于红的像素数（即白缝）。"""
    silhouette = np.clip(layers["blue"] + layers["red"], 0.0, 1.0) > LEVEL
    return int((silhouette & (layers["blue"] < LEVEL) & (layers["red"] < LEVEL)).sum())


def report(rgb: np.ndarray, raw: dict, fields: dict) -> None:
    """打印原稿色彩、连通块与接缝诊断信息。"""
    for name in ("blue", "red"):
        med, mode = color_stats(rgb, raw[name])
        print(f"[{name}] median={med} mode={mode}")
        for area, box in kept_blocks(fields[name]):
            print(f"    block area={area} bbox={box}")
    print(f"[seam] 解混场白缝 {sliver_px(raw)} px（原稿）；重构场 {sliver_px(fields)} px（含外缘软化像素）")


def seam_gap_px(color: np.ndarray) -> int:
    """统计渲染图中夹在蓝红之间的背景像素（真白缝）。"""
    blue = np.linalg.norm(color - rgb_to_float(BLUE_HEX), axis=2) < CLASS_TOL
    red = np.linalg.norm(color - rgb_to_float(RED_HEX), axis=2) < CLASS_TOL
    white = np.linalg.norm(color - np.array(WHITE), axis=2) < CLASS_TOL
    near_b = ndimage.binary_dilation(blue, np.ones((3, 3)))
    near_r = ndimage.binary_dilation(red, np.ones((3, 3)))
    return int((white & near_b & near_r).sum())


def main() -> None:
    """执行矢量化并写出 SVG 与校验图。"""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rgb = load_rgb(SRC)
    raw_blue, raw_red = unmix(rgb, rgb_to_float(BLUE_HEX), rgb_to_float(RED_HEX), np.array(WHITE))
    fields = dict(zip(("blue", "red"), seam_fields(raw_blue, raw_red)))
    report(rgb, {"blue": raw_blue, "red": raw_red}, fields)
    layers = {key: vectorize(value) for key, value in fields.items()}
    layers = {key: [merge_weak_lines(loop) for loop in value] for key, value in layers.items()}
    hull_idx = hull_index(layers["blue"], np.array(STERN_POINT))
    layers["blue"] = tilt_deck(
        layers["blue"], hull_idx, np.array(SAIL_TIP_OLD), DECK_TILT_DEG, DECK_TILT_FALLOFF
    )
    foot = red_foot_point(layers["red"])
    drop = deck_y_at(deck_segments(layers["blue"][hull_idx]), foot[0]) - foot[1] + RED_TOUCH_OVERLAP
    layers["red"] = drop_red_foot(layers["red"], foot[1], drop, RED_DROP_SPAN)
    print(f"[foot] 红色 T 下端下移 {drop:.1f}px 贴到船沿（原缝 {drop - RED_TOUCH_OVERLAP:.1f}px）")
    seams = seam_polylines(fields)
    layers = share_seams(layers, seams, SEAM_MATCH, SEAM_FIT_TOL)
    print(f"[seam] 提取接缝曲线 {len(seams)} 条，蓝红共同对齐到同一条边")
    for key in ("blue", "red"):
        smooth_report(key, layers[key])
    shift = centering_shift(layers, CANVAS)
    SVG_OUT.write_text(build_svg(layers, shift), encoding="utf-8")
    zero = np.zeros(2)
    cov = {key: render_coverage(value, zero) for key, value in layers.items()}
    for key in ("blue", "red"):
        iou, mean, p95, mx = mask_metrics(fields[key], cov[key])
        runs = sum(len(loop) for loop in layers[key])
        segs = sum(len(run) for loop in layers[key] for run in loop)
        print(
            f"[{key}] contours={len(layers[key])} runs={runs} beziers={segs} "
            f"IoU={iou:.4f} dev mean={mean:.3f}px p95={p95:.3f}px max={mx:.3f}px"
        )
    print(f"[layout] centered, shift=({shift[0]:.1f}, {shift[1]:.1f})")
    groups = [(layers["blue"], BLUE_HEX), (layers["red"], RED_HEX)]
    render = render_rgb(groups, zero, "white")
    print(f"[seam] 渲染图中夹在蓝红之间的背景像素 = {seam_gap_px(render)}")
    write_exports(groups, shift)
    print(f"[write] {SVG_OUT}")
    for size in EXPORT_SIZES:
        print(f"[write] {OUT_DIR / f'duivy_logo_{size}.png'}")
    print(f"[write] {OUT_DIR / 'duivy_logo_transparent_1080.png'}")


if __name__ == "__main__":
    main()
