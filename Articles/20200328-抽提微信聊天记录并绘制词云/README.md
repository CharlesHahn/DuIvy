### 1 微信聊天记录数据库破解

#### 1.1 Android微信数据备份

这一步对于某些手机来说可能需要root权限（我的手机很早就没有底裤了，所以……

这里声明我的小破机的型号：魅蓝note3，Android 7.0，Flyme 6.3.0.2A 。

备份微信数据的操作完全在手机里完成，设置->存储与备份，然后在备份列表里选择 软件->微信，这样就会生成一个只含微信数据的备份，并存储在文件根目录的backup文件夹里面；大概的名字应该是com.tencent.mm 。

#### 1.2 找到数据库

把手机backup文件夹里面的com.tencent.mm文件夹一起弄到电脑上。数据库文件EnMicroMsg.db路径大致为：com.tencent.mm->MicroMsg->32位字符的文件夹->EnMicroMsg.db；如果有过多个账号登录微信，那么32位乱码的文件夹可能有多个。数据库文件同目录下还有很多其它数据库，存储着其它的相关信息。

将数据存储在本地所使用的数据库有很多都是sqlite，微信也是，但是微信加密。

#### 1.3 获取数据库密钥

看雪论坛上已经有不少帖子阐述了对于微信的逆向过程。EnMicroMsg.db的加密方式是这样的：首先需要获得手机的IMEI码（一般15位）和微信的uin码（9位），然后将IMEI码和uin码顺序连接，然后用MD5算法加密，小写密文的前七位就是数据库的密码了。

手机IMEI码获取很简单，在手机拨号键盘上输入 `*#06# ` 即可获得，有些手机可能有多个IMEI码，那可以自己每一个都试一下嘛。

uin码的获取稍微有点儿技巧，很多教程可能放弃智取，直接“撞库”了，但是其实，uin码一直我们都有！就在备份出来的 com.tencent.mm 文件夹里，在这个文件夹下搜索`.xml` ，就可以看到很多很多xml文件，找到 auth_info_key_prefs.xml，用文本编辑器打开之后就能看到醒目位置的uin码；如果找不到 auth_info_key_prefs.xml 文件，也没有关系，其它的xml文件里也有uin码；对于微信来说，一个账号对应的uin码应该是唯一的，所以各xml文件里的uin码都一样。打开其它的xml文件，搜索`uin` ，应该就能看到。当然有些xml文件里没有uin码，多找找嘛，反正会有的。

拿到这两串数字之后，利用[MD5加密网站](https://www.sojson.com/encrypt_md5.html)<sup>[1]</sup>（当然也可以其它网站，这个我用着ok）加密，输入`IMEI + uin` （没有+号，俩串字符直接连起来），选取32位小写密文的前7个字符，这就是我们的数据库密钥了。

#### 1.4 打开数据库

利用sqlcipher软件打开sqlite加密数据库。sqlcipher可以通过 clone GitHub 自己编译，但是肯定也有别人分享的exe文件啊，有个公众号叫 `万能搜吧`，关注，回复`导出` 应该就能获得一个[微云分享的link](https://share.weiyun.com/5t5Gywt)<sup>[2]</sup>，这就是windows上可以用的sqlcipher.exe了。我备份了一个，可以[download](http://charles8hahn.pythonanywhere.com/download/sqlcipher.zip)<sup>[3]</sup> 。得到sqlcipher之后，运行它，然后从file->open DB 打开EnMicroMsg.db，提示输入密码，那么我们输入之前得到的7位密钥就可以啦！

#### 1.5 聊天记录导出

sqlcipher提供了两种导出方式，导出成csv或者SQL数据库；导出成csv或许会更方便处理，当然导出成数据库也很好，可以用脚本直接查询数据库抽提数据绘制词云。这里我们导出为csv，打开数据库之后，选择 file->export->Table as csv ，然后我们需要选择导出的table，选择`message` 然后导出。这样就得到csv文件了。

### 2 微信聊天记录词云绘制
#### 2.1 文本预处理

我们得到的csv文件里有一列数据`content` 就是我们的聊天记录了，`content`前面一列应该是我们聊天的对象；聊天记录包含了文本聊天记录，还有一些图片或者红包之类的代码，我们需要对其进行简单的数据清理。csv可以用excel数据透视很方便地进行处理，首先我们可以用数据透视抽提出和某一个人的聊天记录，然后清理掉聊天记录中那些非文本的聊天记录，最后得到一个只包含文本聊天记录的csv文件作为我们绘制词云的原料。
当然这里的处理过程也可以利用脚本自动完成。

关于编码错误：sqlcipher导出的csv可能excel能打开，sublime和python打不开，我们可以切换一下它的编码，利用windows自带的文本编辑器打开csv文件然后选择另存为Unicode编码的csv文件就可以了。

#### 2.2 绘制词云
python绘制词云的的教程有很多，这里就不再赘述。

文本处理：

```python 
# author : charlie 

import jieba     # 中文分词的库
import wordcloud # 词云

# 打开文件
with open("xxxx.csv", 'r', encoding = 'utf-8') as fo:
    content = fo.read()

# 读取聊天记录文本
string = ""
lines = content.strip().strip('\n').split('\n')
error = 0
for line in lines[1:]:
    msg = line.strip().split(',')
    try:
        string += msg[8].strip()
    except:
        error += 1
# msg[8]就是我们的文本了
# 有的时候读取会出错，可能前面数据没搞干净
print(error, len(lines)) 
# 最终error占比还是很少的，我的输出：622 122398

# 分词，处理成wordcloud可用的字符串
words = jieba.cut(string)
deal = ' '.join(words)
```

绘制词云：

```python 
# 初始化wc
wc = wordcloud.WordCloud(
    font_path="simsun.ttf", # 设置中文字体
    background_color = 'white', 
    max_words = 8000, 
    random_state = 20, 
    width = 800,
    height = 600,
    max_font_size = 40)

# 产生词云并保存
pic = wc.generate(deal)
pic.to_file('xxxx.jpg')
```

也可以绘制图片背景的：

```python
import numpy as np
from PIL import Image
from wordcloud import WordCloud, ImageColorGenerator

# 载入背景图片
# 背景图片最好比较简单，轮廓清楚，不要太大，会很慢
coloring = np.array(Image.open("background.jpg"))

# 初始化wc
# 不设max-font-size参数看起来或许更好看？！
wc = WordCloud(
    font_path="simsun.ttf", 
    background_color = 'white', 
    max_words = 8000, 
    random_state = 20, 
    mask = coloring # 添加mask参数
    )
wc.generate(deal)

# 给词云按照图片上色！
image_colors = ImageColorGenerator(coloring)
pic = wc.recolor(color_func=image_colors)

# 保存词云图片
pic.to_file('wordCloudForIvy_final_1.jpg')
```

对于有python的同学来说运行上面的代码自然ok，没有python的同学如果想试试绘制聊天记录的词云，评论区留言，我可以之后尝试打包出exe可执行文件。


### 3 其它参考资料

1. [微信本地DB破解](https://zhuanlan.zhihu.com/p/91257989)<sup>[4]</sup>
2. [PC微信数据库解密](https://bbs.pediy.com/thread-251303.htm)<sup>[5]</sup>
3. [词云绘制](https://amueller.github.io/word_cloud/auto_examples/index.html#example-gallery)<sup>[6]</sup>
4. 推文排版[Markdown Nice](https://www.mdnice.com/)<sup>[7]</sup>
5. 外链索引[wechat format](https://labs.lyric.im/wxformat/)<sup>[8]</sup>

<sub><sup> 看雪论坛牛皮！以前学汇编的时候咋没多去瞧瞧这个论坛orz </sup></sub><br/>
<sub><sup> 本文封面是我做的另一个资料的词云，聊天自然不会这么多学科交叉啊，捂脸~ </sup></sub><br/>

### 4 References

```
[1] MD5加密网站: 
    https://www.sojson.com/encrypt_md5.html

[2] 微云分享的link: 
    https://share.weiyun.com/5t5Gywt

[3] download: 
    http://charles8hahn.pythonanywhere.com/download/sqlcipher.zip

[4] 微信本地DB破解: 
    https://zhuanlan.zhihu.com/p/91257989

[5] PC微信数据库解密: 
    https://bbs.pediy.com/thread-251303.htm

[6] 词云绘制: 
    https://amueller.github.io/word_cloud/auto_examples/index.html

[7] Markdown Nice: 
    https://www.mdnice.com/

[8] wechat format: 
    https://labs.lyric.im/wxformat/
```








