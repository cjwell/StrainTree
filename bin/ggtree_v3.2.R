
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
suppressMessages(library(treeio))
suppressMessages(library(ape))
suppressMessages(library(phytools))
suppressMessages(library(getopt))

arg <- matrix(
  c("input", "i", 1, "character", "输入newick格式的进化树文件",
  "output", "o", 2, "character", "输出文件名称,无需扩展名[省略时将使用输入文件名]",
  "layout", "l", 2, "character", "进化树布局[默认:rectangular], 可选(rectangular, slanted, ellipse, roundrect\n\t\t\t       circular, fan, radial, ape, daylight,)",
  "label_branch", "b", 2, "character", "排列标签(tips)或忽略枝长[a:排列标签, b:忽略枝长,默认(省略):不排列标签不忽略枝长]",
  "font_size", "f", 2, "numeric", "指定标签(tips)文字的字体大小[默认值:2.3](过大的字体可能会导致tips名称超出图形范围)\n\t\t\t    (当tips名称超出图形范围时也可以通过缩小字体以达到完全显示)",
  "line_size", "s", 2, "numeric", "指定进化树的线宽[默认值:0.4]",
  "point_size", "p", 2, "numeric", "指定bootstrap value 节点的点大小[默认值:2]",
  "node_number", "w", 2, "character", "是否显示节点位置(蓝色数字)[a:显示全部节点,i:仅显示内部节点]",
  "angle", "a", 2, "numeric", "指定圆形树形的开口角度(仅当树形布局 --layout | -l 为fan时生效)",
  "bootstrap_value","v", 2, "character", "是否在内部节点处显示bootstrap的值(紫色数值)[默认:不显示,输入数值将会调整数值的水平位置默认为1.2]",
  "outgroup_s", "g", 2, "character", "输入需要当作外枝的tips名称(关键词),用以reroot进化树(仅当外枝为一个tips时)",
  "outgroup_m", "j", 2, "character", "提供某一内部节点后包含的tips名称列表,以在该内部节点处reroot进化树(每行一个tips)",
  "reroot_node", "t", 2, "integer", "选择一个节点用于reroot进化树(可以使用 -w 参数先显示节点位置,以确定要reroot的节点)",
  "rename_tips", "m", 2, "character", "输入要进行替换的tips名称列表(列表需包含原名和替换名两列并用制表符分隔)",
  "name_length", "n", 2, "integer", "保留tips名称的长度(字符数)",
  "root_length", "r", 2, "numeric", "在根结点处添加枝长,该值将乘以平均枝长(通常作用于圆形树形以增大圆形的内部空间)",
  "show_bootstrap", "z", 2, "character", "是否显示bootstrap[默认:T 显示, F:不显示]",
  "width", "d", 2, "numeric", "等比例调整输出pdf的宽度(默认值:1)",
  "height", "e", 2, "numeric", "等比例调整输出pdf的高度(默认值:1)",
  "help", "h", 0, "logical", "输出帮助信息!"
  ), byrow = TRUE, ncol = 5)
opt <- getopt(spec = arg)
input <- opt$input
output <- opt$output
line_size <- opt$line_size
font_size <- opt$font_size
layout<- opt$layout
label_branch <- opt$label_branch
point_size <- opt$point_size
root_length <- opt$root_length
node_number <- opt$node_number
show_bootstrap <- opt$show_bootstrap
outgroup_s <- opt$outgroup_s
outgroup_m <- opt$outgroup_m
rename_tips <- opt$rename_tips
reroot_node <- opt$reroot_node
bootstrap_value <- opt$bootstrap_value
angle <- opt$angle
name_length <- opt$name_length
hjpdf <- opt$height
wjpdf <- opt$width

if (!is.null(opt$help) || is.null(opt$input)) {
  cat("\n使用方法:", 
  paste("\n", getopt(spec = arg, usage = TRUE), "\n"))
  quit()
}

#初始化值
if (is.null(point_size)) {
  point_size <- 2
}
if (is.null(output)) {
output <- input
}
if (is.null(angle)) {
  angle <- 0
}
if (is.null(line_size)) {
  line_size <- 0.4
}

if (is.null(font_size)) {
  font_size <- 2
}

if (is.null(layout)) {
  layout <- "rectangular"
}
if(is.null(hjpdf)) {
  hjpdf <- 1
}
if( is.null(wjpdf)) {
  wjpdf <- 1
}
# 忽略枝长或者对齐标签
if (is.null(label_branch)) {
  align_tips <- FALSE
  branch_length <- "branch.length"
}else if (label_branch == "a") {
  align_tips <- TRUE
  branch_length <- "branch.length"
}else if (label_branch == "b") {
  align_tips <- FALSE
  branch_length <- "none"
}

#读入文件
treedata <- treeio::read.newick(input, node.label = "support")
#将名称中的-替换成下划线
treedata@phylo$tip.label <- gsub("-",'_',treedata@phylo$tip.label)


df<- fortify(treedata)
ntip <- length(treedata@phylo$tip.label)
legendsize <- (8 / 110) * ntip + 48 / 11

#根据输入的需要reroot的tips的名称推断其根节点
##当提供一个文件列表时
if (!is.null(outgroup_m)){
  outgroupm <- read.table(outgroup_m, header = FALSE,sep = "\t")
  outgroupm[,1] <- gsub("-", "_", outgroupm[,1])
  outgroupm[,1] <- gsub(" ", "", outgroupm[,1])
  a<- vector()
  for (i in 1:length(outgroupm[,1])){
    a[i] <- as.integer(df[grep(outgroupm[i,1],df$label), 2])
  }
  treedata_phylo_groupm <- ape::as.phylo(treedata)
  b <- vector()
  for (t in 1:length(a)){
    b[t] <- as.integer(phytools:::getParent(treedata_phylo_groupm, node = a[t]))
  }
  reroot_node <- min(b)
}

#当名称为只有一个时
if (!is.null(outgroup_s)) {
outgroup_s <- gsub("-", "_",outgroup_s)
reroot_node <- as.integer(df[grep(outgroup_s,df$label, fixed =TRUE), 2])
}

#tips名称替换
if (!is.null(rename_tips)) {
  
    name <- read.table(rename_tips, header = FALSE,sep = "\t")
	name[1] <- gsub("-", "_",name[,1],fixed = TRUE )
    colnames(name) <- c("ori", "new")
    for (i in name$ori) {

      treedata@phylo$tip.label[grep (i, treedata@phylo$tip.label)] <- name[which(name$ori == i), 2]
    }
 
}

new_node <- (getNodeNum(treedata) + 1)
rebranch <- treedata@phylo$edge.length[which(treedata@phylo$edge[,2] == reroot_node)]/2

#reroot进化树
if (is.null(reroot_node)) {
  q<- 1
} else if (reroot_node <= ntip) {	#在tips节点处reroot进化树
  rootnode <- rootnode(treedata)
  treedata_phylo <- ape::as.phylo(treedata)
  parent_node <- phytools::getParent(treedata_phylo, node = reroot_node)
  treedata_phylo_reroot<- phytools::reroot(treedata_phylo, node.number = reroot_node, position = rebranch )
  treedata@phylo <- treedata_phylo_reroot
  treedata@data$node[which(treedata@data$node == parent_node)] <- NA
  treedata@data$node[which(treedata@data$node< parent_node & treedata@data$node >rootnode)] <- treedata@data$node[which(treedata@data$node < parent_node & treedata@data$node > rootnode)]+1
} else if (reroot_node > ntip) {	#在内部节点处reroot进化树
  newnode_support <- data.frame(new_node,treedata@data$support[which(treedata@data$node == reroot_node)])
  colnames(newnode_support) <- c('node','support')
  treedata@data$node[which(treedata@data$node==reroot_node)]<- new_node
  treedata<- ape::root(treedata,node=reroot_node,resolve.root = TRUE,)
  treedata@phylo$edge.length[which(treedata@phylo$edge[,2]==reroot_node)] <- rebranch
  treedata@phylo$edge.length[which(treedata@phylo$edge[,2]==new_node)] <- rebranch
  treedata@data <- rbind(treedata@data,newnode_support)
}

#统一bootstrap value的值
if(max(treedata@data$support,na.rm = TRUE) <= 1) {
  treedata@data$support<- treedata@data$support * 100
}

#ggtree绘图组件
tree <- ggtree(treedata, size = line_size, layout = layout, branch.length = branch_length, open.angle = angle)
tiplabel <- geom_tiplab(size = font_size,align = align_tips)
scale <- geom_treescale(x = 0, y = -2, linesize = line_size, fontsize = font_size + 0.5, color = "black")
nodepoint <- geom_nodepoint(aes(subset = (!isTip & !is.na(support)), 
                                  fill = cut(support, c(0, 70, 90, 100))),
                              size = point_size, shape = 21, na.rm = TRUE)
fill_manual<- scale_fill_manual(values = c("black", "grey", "white"), 
                  guide ='legend', name = 'Bootstrap Percentage(BP)', 
                  breaks = c('(90,100]', '(70,90]', '(0,70]'), 
                  labels = expression(BP>=90,70 <= BP * " < 90", BP < 70))
tree_theme <- theme_tree(legend.position = "bottom",
                         plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
                         plot.title = element_text(size = legendsize+4, hjust = 0.5),
						 legend.text = element_text(size = legendsize), 
                         legend.title = element_text(size = legendsize))
coord <- coord_cartesian(clip = "off")
title <- ggtitle(input)

#限制名称长度
if (!is.null(name_length)) {
  tree$data$label[which(tree$data$isTip == TRUE)] <- substr(tree$data$label[which(tree$data$isTip == TRUE)], 0, name_length)
}

#是否显示node的位置
if (is.null(node_number)) {
  node_position <- xlab("")
} else if (node_number == "i") {
  node_position <- geom_text2(aes(label = node, subset = !isTip), color = "#2986CC", size = font_size, hjust = -0.3)
} else if (node_number == "a") {
  node_position <- geom_text2(aes(label = node), color = "#2986CC", size = font_size, hjust= -0.3)
}

#是否显示bootstrap的值
if (!is.null(bootstrap_value)) {
  if(bootstrap_value == 'T' ){
  bootstrap_value <- 1.2
    show_support <- geom_text2(aes(label=support), color = "#C90076",size = font_size,hjust = bootstrap_value,vjust = -0.4)
  } else if (bootstrap_value == 'F') {
    show_support <- xlab('')
  } else {
  bootstrap_value <- as.numeric(bootstrap_value)
  show_support <- geom_text2(aes(label=support), color = "#C90076",size = font_size,hjust = bootstrap_value,vjust = -0.4)
  }
} else {
  show_support <- xlab('')
}

#提取树信息用于限制图形大小
maxn<- max(nchar(tree$data$label[which(tree$data$isTip == TRUE)]))
maxh <- length(tree$data$label[which(tree$data$isTip == TRUE)]) 
maxw <- max(tree$data$branch) 

#设置根的枝长
if (!is.null(reroot_node)) {

  if(!is.null(root_length)) {
    root_length <- root_length
  } else {
    root_length <- 1
  }
  rootedge <- geom_rootedge(rebranch * root_length, size = line_size)
} else {
  if(!is.null(root_length)) {
    root_length <- root_length
  } else {
    root_length <- 0
  }
  rootedge <- geom_rootedge(mean(tree$data$branch.length) * root_length, size = line_size)
}

#设置是否以点的形式显示bootstrap
if(is.null(show_bootstrap) || show_bootstrap=='T') {
  p<- 1
}else if(show_bootstrap == "F") {
  nodepoint <- xlab('')
  fill_manual<- xlab('')
}

## 输出图形类型判断
if (branch_length != "none") {

  if (layout == "rectangular" | layout == "ellipse" | layout == "slanted" | layout == "roundrect") {

	pdfh <- (0.1091 * maxh) + 3
    pdfw <- pdfh * (0.994^maxh + 0.25)
	pdfh <- pdfh *hjpdf
	pdfw <- pdfw *wjpdf
	hjust_title<- (((maxn/7.5)/2.54)/(pdfw-0.5-((maxn/7.5)/2.54)))*0.5 + 0.5
    tree_theme <- theme_tree(plot.margin = margin(t = 0.5, r = maxn/7.5, b = 0.5, l = 0.5, unit = "cm"),
                             legend.position = "bottom", legend.direction = "horizontal",
							 legend.key.size = unit(0.1, "cm"),
                             legend.text = element_text(size = legendsize), 
                             legend.title = element_text(size = legendsize), 
                             plot.title = element_text(size = legendsize + 4, hjust = hjust_title))

    p <- tree + tiplabel + nodepoint + tree_theme + fill_manual + node_position + title + 
	    rootedge + show_support+coord + scale

  }else if (layout == "circular" | layout == "fan" | layout == "inward_circular" | layout == "radial") {
    pdfh <- maxh / 25 + 8
    pdfw <- maxh / 25 + 8
	pdfh <- pdfh *hjpdf
	pdfw <- pdfw *wjpdf
    expandh <- hexpand(maxn / maxh, direction = 1)
    expandv <- vexpand(maxn / maxh, direction = -1)
    
    p <- tree + tiplabel + nodepoint + tree_theme + fill_manual + expandh + expandv + title + 
       node_position + scale + show_support +rootedge
	  	ggsave(p , filename = paste(output, ".pdf", sep = ""),width = pdfw, height = pdfh, limitsize = FALSE)
  }else if (layout == "ape" | layout == "daylight" | layout == "equal_angle") {
    pdfh <- maxh / 25 + 8
    pdfw <- maxh / 25 + 8
	pdfh <- pdfh *hjpdf
	pdfw <- pdfw *wjpdf
    expandh <- ggexpand(maxn / 100, direction = -1, side = "hv")
    expandv <- ggexpand(maxn / 100, direction = 1, side = "hv")
    p <- tree + tiplabel + nodepoint + fill_manual + tree_theme + coord + title + expandh +
      expandv + node_position + rootedge +show_support
  }else if (layout == "dendrogram") {
    pdfw <- (0.1091 * maxh) + 2.691
    pdfh <- pdfw * (0.994^maxh + 0.2)
	pdfh <- pdfh *hjpdf
	pdfw <- pdfw *wjpdf
	theme_den <- theme_dendrogram(plot.margin=margin(1,1,maxn/7.5,1,unit = "cm"),
	                              plot.title = element_text(size = legendsize+4, hjust = 0.5),
								  legend.text = element_text(size = legendsize),
                                  legend.title = element_text(size = legendsize))
    tiplabel <- geom_tiplab(size = font_size, align = align_tips, color = "black", vjust = 0.5, hjust = 1, angle = 90)
    p <- tree + tiplabel + nodepoint + fill_manual + tree_theme + title + node_position + show_support+theme_den
  }
}else if (branch_length == "none") {
  
  if (layout == "rectangular" | layout == "ellipse" | layout == "slanted" | layout == "roundrect") {
    pdfh <- (0.1091 * maxh) + 3
    pdfw <- pdfh * (0.994^maxh + 0.2)
	pdfh <- pdfh *hjpdf
	pdfw <- pdfw *wjpdf
	hjust_title<- (((maxn/7.5)/2.54)/(pdfw-0.5-((maxn/7.5)/2.54)))*0.5 + 0.5
    tree_theme <- theme_tree(plot.margin = margin(t = 1, r = maxn / 7, b = 1, l = 1, unit = "cm"),
                         legend.position = "bottom", 
                         legend.text = element_text(size = legendsize),
                         legend.title = element_text(size = legendsize), 
                         plot.title = element_text(size = legendsize + 4, hjust = hjust_title), 
                         legend.key.size = unit(0.1, "cm"))
    
    p <- tree + tiplabel + scale + nodepoint + tree_theme + fill_manual+ 
	  coord + title + node_position + show_support
  }else if (layout == "circular" | layout == "fan" | layout == "inward_circular" | layout == "radial") {
    pdfh <- (maxh / 25 + 8)
    pdfw <- (maxh / 25 + 8)
	pdfh <- pdfh *hjpdf
	pdfw <- pdfw *wjpdf
    expandh <- hexpand(maxn / maxh, direction = 1)
    expandv <- vexpand(maxn / maxh, direction = -1)
    
    p <- tree + tiplabel + nodepoint + tree_theme + fill_manual + expandh + expandv + 
	  title + node_position + show_support
  }else if (layout == "ape" | layout == "daylight" | layout == "equal_angle") {
    pdfh <- (maxh / 25 + 8)
    pdfw <- (maxh / 25 + 8)
	pdfh <- pdfh *hjpdf
	pdfw <- pdfw *wjpdf
    expandh <- ggexpand(maxn / 100, direction = -1, side = "hv")
    expandv <- ggexpand(maxn / 100, direction = 1, side = "hv")
    p <- tree + tiplabel + tree_theme + nodepoint + fill_manual+ coord + title + expandh +
      expandv + node_position + show_support
  }else if (layout == "dendrogram") {
    pdfw <- (0.1091 * maxh) + 2.691
    pdfh <- pdfw * (0.994^maxh + 0.2)
	pdfh <- pdfh *hjpdf
	pdfw <- pdfw *wjpdf
	theme_den <- theme_dendrogram(plot.margin=margin(1,1,maxn/7.5,1,unit = "cm"),
	                              plot.title = element_text(size = legendsize+4, hjust = 0.5),
								  legend.text = element_text(size = legendsize),
                                  legend.title = element_text(size = legendsize))
    tiplabel <- geom_tiplab(size = 2, align = align_tips, color = "black", vjust = 0.5, hjust = 1,angle = 90)
    p <- tree + tiplabel + nodepoint + fill_manual + tree_theme + 
	title + node_position + show_support +theme_den 
  }
}
	

ggsave(p , filename = paste(output, ".pdf", sep = ""), width = pdfw, height = pdfh, limitsize = FALSE)