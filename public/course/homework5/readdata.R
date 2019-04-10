# 下面命令用于查看目前工作目录，然后把数据
# data.RData放进该工作目录
getwd() 
# 导入数据，如果该数据没有放入工作目录，
# 则需要用绝对路径，或者用setwd()更改工作目录
# 至数据集所在文件夹
load("data.RData") 
# 便于直接使用该数据集的两个变量height和sex用attach命令
attach(heights) 