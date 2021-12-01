# docsify and its application

## docsify 简介

[docsify](https://docsify.js.org/) 官方介绍是"a magical documentation site generator"，简言之就是可以从markdown文档直接生成文档网站，且不会生成静态的html文件，而是在运行时对md文档进行转换。它非常的轻量，同时提供了多种插件来进行功能拓展，还有智能的全文搜索，美滋滋！

具体情况和tutorial可以参考docsify的官网和github。

安装方法的话，可以先安装nodejs和npm，再按照官网的方法进行全局安装。如果没有npm，也可以按照官网的快速开始的手动初始化方法进行初始化。

初始化之后会有一个index.html，通常还需要如下文件：

- _coverpage.md 自定义封面用的
- _sidebar.md 自定义目录栏用的
- README.md 自定义进去的第一页
- 各种插件的js文件（可以直接在index.html里面引用，也可以下载到本地）

所有这些如何生成，官网都有细致的描述。

然后把自己的markdown文档放在同一目录下或者更深一级目录，在_sidebar.md里面写好链接就可以啦。

## 结合 github pages 建站

我们生成了这样的一个文件夹，包含docsify需要的几个文件和我们的文档，可以用git来管理，然后推送到github或者gitee，然后利用github或者gitee的pages功能来生成互联网上可访问的公开的文档网站。具体操作可以参照 https://opensource.com/article/20/7/docsify-github-pages 或者自行STFW。



## 结合云主机建站

可以参考网上的众多教程。

或者可以偷懒，直接把docsify的这个文件夹copy到云主机上，在云主机上安装npm和docsify然后`docsify serve your_docs/`，或者在文件夹路径下`python -m http.server 3000`，如果你的云主机有公网ip，这就足够了，访问的话就访问相应的公网ip加端口号，如 xxx.xxx.xx.xx:3000

要是在自己电脑上预览了，那么同一局域网下的其它设备也是可以通过ip:port的方式访问的。

## 结合 pythonanywhere 建站

当然还可以结合pythonanywhere和flask（或者Django）来建站。pythonanywhere一个邮箱就可以免费注册，并且做好的网站可以在互联网上访问；相比github pages，不方便的就是需要自己手写下网站后端代码，对于docsify这种文档网站来说，十几行代码就足够了；还有一个不方便的就是可能不能像github那样通过git及时同步，需要手动操作，就很烦，但这个问题也可以通过写一个同步控制脚本完成（模拟登录pythonanywhere然后把文档上传到相应文件夹就行了，没啥难度）。



## something else

我用git管理自己的文档，写作和记笔记主要是neovim，docsify主要方便自己展示和搜索信息来着；当然，也推荐[the_silver_searcher](https://github.com/ggreer/the_silver_searcher)，一个grep的替代，windows10也有的用，用来做文档内容搜索很不错；linux下的pdf搜索可以用pdfgrep，也好用，不过还没发现win10有编译好的，如有请推荐。

上文很简略，本想详尽写写，但没必要，有需要的人自然会去STFW或者RTFM；docsify和github都很简单，比较困难的可能就pythonanywhere的后端代码和目录结构，自己研究研究吧，很有裨益的一件事情。又或者，哪天再来重写一份非常详细的。

18年夏就念叨着说要学习julia，前段时间终于开始认真把自己的一些代码用julia重写了。总结起来，现今的julia已经比两年前好用多了，如果有兴趣，其实可以考虑入手。当然，先要搞清楚自己的目的。python在脚本、爬虫等等领域还是非常有优势的，在机器学习等领域，python也有大量丰富的第三方库，Julia也不及。倒是数值计算领域，如果嫌弃python慢，可以考虑转向Julia。要小心，个人感觉Julia并不是很好学，太多的特性，号称学习了matlab、python、R、C等语言的优点，其实语法有点四不像；入手还是有一些难度（相较于python），也并不是小白就可以写出高性能Julia代码的（也别想用Julia写各种小脚本，光是JIT的预热启动就够让人泄气的了）。

上个月也把键盘布局从qwerty迁移到了colemak，倒是成功减轻了我长期使用键盘的手部疼痛，但是到了别人的键盘上就秒变白痴，有利有弊吧。

最近科研没思路，烦躁，脾气不好。

祝好~

