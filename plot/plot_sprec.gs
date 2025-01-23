"reinit"
"set display color white"
"c"

"set lwid 13 6"
"set lwid 14 12"
"set lwid 15 4"
"set annot 1 15"
"set strsiz 0.19"
"set xlopts 1 15 0.25"
"set ylopts 1 15 0.25"
"set clopts 1 15 0.25"
"set rgb 200 100 100 100 220"
"set grid on 3 200 6"

"set font 11 file /home/der0318/.grads/Helvetica.ttf"
"set font 11"

path="/data/der0318/taiwanVVM/"

nrun=40
run.1='tpe20050702nor'
run.2='tpe20050712nor'
run.3='tpe20050723nor'
run.4='tpe20050826nor'
run.5='tpe20060508nor'
run.6='tpe20060623nor'
run.7='tpe20060702nor'
run.8='tpe20060718nor'
run.9='tpe20060721nor'
run.10='tpe20070830nor'
run.11='tpe20080715nor'
run.12='tpe20090707nor'
run.13='tpe20090817nor'
run.14='tpe20090827nor'
run.15='tpe20100629nor'
run.16='tpe20100630nor'
run.17='tpe20100802nor'
run.18='tpe20100803nor'
run.19='tpe20100912nor'
run.20='tpe20110615nor'
run.21='tpe20110616nor'
run.22='tpe20110702nor'
run.23='tpe20110723nor'
run.24='tpe20110802nor'
run.25='tpe20110816nor'
run.26='tpe20110821nor'
run.27='tpe20120715nor'
run.28='tpe20120819nor'
run.29='tpe20130630nor'
run.30='tpe20130703nor'
run.31='tpe20130705nor'
run.32='tpe20130723nor'
run.33='tpe20130807nor'
run.34='tpe20130825nor'
run.35='tpe20140525nor'
run.36='tpe20140703nor'
run.37='tpe20140711nor'
run.38='tpe20140714nor'
run.39='tpe20140825nor'
run.40='tpe20150613nor'


i=1
while(i<=nrun)

"open "path"/"run.i"/gs_ctl_files/surface.ctl"
"open "path"/"run.i"/gs_ctl_files/topo.ctl"

"set mpdraw off"

"set grads off"
"set parea 1.2 10 1.2 7.5"
"set ylint 1"
"set xlint 1"

"color -levs 1 2 6 10 15 20 30 40 50 70 90 110 130 150 200 300 -kind lightgray->(160,255,195)->(0,205,255)->(0,150,255)->(0,105,250)->(50,150,0)->(50,255,0)->(255,255,0)->(255,200,0)->(255,150,0)->(255,0,0)->(200,0,0)->(160,0,0)->(150,0,155)->(200,0,210)->(255,0,245)->(255,100,255)->(255,200,255)"
"set gxout grfill"
"d sum(sprec,t=1,t=49)*600." 
*"d sum(sprec,t=1,t=145)*600." 

"set line 1 1 15"
"xcbar 8.65 8.8 1.7 7 -edge triangle -line on -fw 0.25 -fh 0.25"

"set string 1 c 15 0"
"set strsiz 0.3"
"draw string 5.6 7.8 "run.i""

"set gxout contour"
"set cthick 15"
"set ccolor 1"
"set clevs 0.05"
"set clab off"
"d topo.2"

"set cthick 10"
"set rgb 183 80 80 80 200"
"set ccolor 183"
"set clevs 0.5"
"set clab off"
"d topo.2"

"printim ./figure/morning08_pr_"run.i".png x2048 y1536"
*"printim ./figure/day_pr_"run.i".png x2048 y1536"
"c"

"close 2"
"close 1"

i=i+1
endwhile



