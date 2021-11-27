# EReg(SPSS) 2021.11.27
## Extending regression analysis

EReg是我为了扩展SPSS在回归方面的功能而开发的插件，采用Apache 2.0开源许可证。坏消息是，由于单人开发，个人水平和精力都有限，加上缺乏写代码的经验，EReg的语法注释较少，存在较多可以改进的地方（当前少有人使用SPSS代码编写程序，短时间内可能不会对代码进行优化）。好消息是，基本功能都有，操作简单且有详细的用户手册（本手册）。

EReg现包含绘图（visualization,即可视化）、多元回归（multiple regression）、二次效应（quadratic effect）、调节效应（moderation effect）、样条回归（spline regression）、岭回归（ridge regression）和主成分回归（principal component regression）等七个模块。原本计划开发偏最小二乘模块，但该方法一般出现在多个因变量的背景下，与EReg中现有模块相性较差，就被放弃了。尽管有些模块的使用频率可能不高，但还是本着“可以不用，但不能没有”的想法开发出来了，毕竟使用率低不一定是方法不好，也可能是方法尚未普及。希望这些功能对用户有所帮助。除了EReg的分析功能，用户手册内提供了大量参考文献以及其它软件的代码，希望这些工作也能发挥其作用。EReg的更新信息和扩展视频见微信公众号（邱宗满，ID：qiuzongman2020）与Bilibili（邱宗满，ID：423767625）：

![微信公众号](https://github.com/zongmanqiu/EReg/blob/main/%E5%9B%BE/%E5%BE%AE%E4%BF%A1%E5%85%AC%E4%BC%97%E5%8F%B7.jpg)
![Bilibili](https://github.com/zongmanqiu/EReg/blob/main/%E5%9B%BE/Bilibili.jpg)

EReg是给SPSS编写的扩展插件，但我很希望在什么时候有国内开发者也做一个与SPSS类似的窗口化操作软件，能开源，能集成Python和R，能使用命令行分析，对中文使用者也足够友好。如果有这样的软件，那一定很令人激动。
