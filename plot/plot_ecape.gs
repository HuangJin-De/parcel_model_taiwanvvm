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
"open ../gs_ctl/cape_"run.n".ctl"
"open "path"/"run.n"/gs_ctl_files/topo.ctl"

"set mpdraw off"

"set grads off"
"set parea 1.2 5.4 3.2 7.4"
"set xlint 1"
"set ylint 1"

"color 0 600 50 -kind white->p31->p32->p33->p34->p35->p37->p38->p39->p310->p311"
"set gxout grfill"
"d ecape2(t=1)"

"set gxout contour"
"set cthick 15"
"set ccolor 1"
"set clevs 0.05"
"set clab off"
"d height.2"

"set cthick 10"
"set rgb 183 80 80 80 200"
"set ccolor 183"
"set clevs 0.5"
"set clab off"
"d height.2"




"set grads off"
"set parea 6.0 10.2 3.2 7.4"
"set xlint 1"
"set ylint 1"

"color 0 600 50 -kind white->p31->p32->p33->p34->p35->p37->p38->p39->p310->p311"
"set gxout grfill"
"d ecape2(t=49)"

"set gxout contour"
"set cthick 15"
"set ccolor 1"
"set clevs 0.05"
"set clab off"
"d height.2"

"set cthick 10"
"set rgb 183 80 80 80 200"
"set ccolor 183"
"set clevs 0.5"
"set clab off"
"d height.2"



"set line 1 1 15"
"xcbar 1.2 10.2 2.3 2.5 -line on -fs 4 -fw 0.2 -fh 0.2"

"set string 1 c 15 0"
"set strsiz 0.3"
"draw string 5.65 7.9 "run.n""

"printim ./figure/decape_"run.n".png x2048 y1536"
"c"

"close 2"
"close 1"
n=n+1
endwhile

"quit"
