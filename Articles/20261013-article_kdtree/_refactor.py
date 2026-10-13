# -*- coding: utf-8 -*-
"""将推文 md 的链接按 GB/T 7714 顺序编码制转为参考文献引用。"""
import re

path = r"C:\Users\hhhhh\Desktop\DuIvyInteraction\article_kdtree\KDTree距离计算优化_推文.md"
md = open(path, encoding="utf-8").read()

# 编号表（首次出现顺序）
refs = [
    ("https://en.wikipedia.org/wiki/K-d_tree", "K-d tree. Wikipedia. https://en.wikipedia.org/wiki/K-d_tree"),
    ("https://cstheory.stackexchange.com/questions/46235/time-complexity-for-nearest-neighbor-searches-in-kd-trees",
     "Time complexity for NN searches in kd-trees. StackExchange (cstheory). https://cstheory.stackexchange.com/questions/46235/time-complexity-for-nearest-neighbor-searches-in-kd-trees"),
    ("https://www.cs.cornell.edu/courses/cs4780/2023sp/lectures/lecturenote16.html",
     "KD Trees. Cornell CS4780 (2023 fall). https://www.cs.cornell.edu/courses/cs4780/2023sp/lectures/lecturenote16.html"),
    ("https://faculty.cc.gatech.edu/~isbell/reading/papers/cover-tree-icml.pdf",
     "Beygelzimer A, Kakade S, Langford J. Cover Trees for Nearest Neighbor. ICML, 2006. https://faculty.cc.gatech.edu/~isbell/reading/papers/cover-tree-icml.pdf"),
    ("https://stackoverflow.com/questions/35930246/range-search-complexity-of-r-tree-and-r-tree",
     "R-tree 范围查询复杂度. StackOverflow. https://stackoverflow.com/questions/35930246/range-search-complexity-of-r-tree-and-r-tree"),
    ("https://arxiv.org/html/2102.08942v1",
     "A Survey on Locality Sensitive Hashing. arXiv:2102.08942. https://arxiv.org/html/2102.08942v1"),
    ("https://www.mongodb.com/resources/basics/hierarchical-navigable-small-world",
     "What is a Hierarchical Navigable Small World (HNSW). MongoDB. https://www.mongodb.com/resources/basics/hierarchical-navigable-small-world"),
    ("https://icml.cc/Conferences/2009/papers/360.pdf",
     "Radovanović M, Nanopoulos A, Ivanović M. Nearest Neighbors in High-Dimensional Data: The Emergence and Influence of Hubs. ICML, 2009. https://icml.cc/Conferences/2009/papers/360.pdf"),
    ("https://pyimagesearch.com/2025/01/27/approximate-nearest-neighbor-with-locality-sensitive-hashing-lsh/",
     "Implementing Approximate Nearest Neighbor Search with LSH. PyImageSearch. https://pyimagesearch.com/2025/01/27/approximate-nearest-neighbor-with-locality-sensitive-hashing-lsh/"),
    ("https://www.pinecone.io/learn/series/faiss/hnsw/",
     "HNSW. Pinecone Learning. https://www.pinecone.io/learn/series/faiss/hnsw/"),
]
# url → 编号
url2num = {url: i + 1 for i, (url, _) in enumerate(refs)}
# 也映射 2018fa 版 Cornell → 归 2023sp [3]（同源）
url2num["https://www.cs.cornell.edu/courses/cs4780/2018fa/lectures/lecturenote16.html"] = 3

# 1. 正文链接 → [n]
def repl(m):
    text, url = m.group(1), m.group(2)
    num = url2num.get(url.lstrip())
    if num:
        return f"{text}[{num}]"
    # 找不到（未预期链接）保留原文
    return m.group(0)

# 只处理正文（参考资料段之前）
ref_marker = "**参考资料**"
body, refs_part = md.split(ref_marker, 1)
body = re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)', repl, body)
md_new = body + ref_marker + refs_part

# 2. 文末参考资料段：替换为 10 条
ref_block = "**参考资料**\n\n" + "\n".join(f"{i+1}. {desc}" for i, (_, desc) in enumerate(refs)) + "\n"
# 找到原参考资料段到 --- 之间替换
md_new = re.sub(r'\*\*参考资料\*\*.*?\n---', ref_block + "\n---", md_new, flags=re.S)

open(path, "w", encoding="utf-8").write(md_new)
print("=== 完成 ===")

# 验证
md2 = open(path, encoding="utf-8").read()
links_left = re.findall(r'\[([^\]]+)\]\(https?://', md2)
print("剩余链接（应 0）:", len(links_left))
print("引用 [n] 数:", len(re.findall(r'\[\d+\]', md2)))
refs_final = md2[md2.find("**参考资料**"):]
print("\n=== 新文末参考资料 ===")
print(refs_final[:800])