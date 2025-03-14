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

path="/data/der0318/parcel_model_taiwanvvm/taiwanVVM/"


lev.1=4000
lev.2=5000
lev.3=6000
lev.4=7000

parea.1="2.0 5.0 4.6 7.6"
parea.2="5.5 8.5 4.6 7.6"
parea.3="2.0 5.0 1.4 4.4"
parea.4="5.5 8.5 1.4 4.4"

nrun=sys("ls "path" | wc -l")
nrun=sublin(nrun,1)

namelist=sys("ls "path)
n=1
while (n<=nrun)
run.n=subwrd(namelist,n)
n=n+1
endwhile


n=1
while(n<=nrun)
"open "path"/"run.n"/gs_ctl_files/dynamic.ctl"
"open "path"/"run.n"/gs_ctl_files/topo.ctl"

"set t 49"
"set mpdraw off"

i=1
while(i<=4)

"set lev "lev.i""
"set parea "parea.i""

"set grads off"
"set xlab off"
"set ylab off"

"color 0 10 1 -kind white->p31->p32->p33->p34->p36->p37->p38->p39->p310->p311"
"set gxout grfill"
"d mag(u,v)"

"set line 1 1 15"
"xcbar 8.7 8.9 1.4 7.6 -line on -fs 2 -fw 0.15 -fh 0.15"

"set gxout vector"
"set cthick 15"
"set ccolor 1"
"set arrowhead 0.06"
"set arrscl 0.5 10"
"set arrlab off"
"d skip(u,128,128);v"


"set gxout contour"
"set cthick 15"
"set ccolor 1"
"set clevs 0.05"
"set clab off"
"d topo.2(z=1,t=1)"

"set cthick 10"
"set rgb 183 80 80 80 200"
"set ccolor 183"
"set clevs 0.5"
"set clab off"
"d topo.2(z=1,t=1)"

i=i+1
endwhile

"set string 1 c 15 0"
"set strsiz 0.3"
"draw string 5.25 7.9 "run.n""

"printim ./figure/wind_"run.n".png x2048 y1536"
"c"

"close 2"
"close 1"
n=n+1
endwhile


