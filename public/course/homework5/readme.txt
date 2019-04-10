getwd() #该命令用于查看目前工作目录，然后把数据data.RData放进该工作目录
load("data.RData") #导入数据，如果该数据没有放入工作目录，则需要用绝对路径，或者用setwd()更改工作目录至数据集所在文件夹
attach(heights) #便于直接使用该数据集的两个变量height和sex
