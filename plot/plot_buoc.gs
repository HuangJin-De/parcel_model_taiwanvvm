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
"open ../gs_ctl/buoc_"run.n".ctl"

"set lev 0 12000"

"set grads off"
"set parea 1.3 10 1.3 7.5"
"set ylabs 0|2|4|6|8|10|12"
"set vrange -0.2 0.2"

"set ccolor 2"
"set cthick 15"
"set cmark 0"
"d b(t=1)"
"set ccolor 11"
"set cthick 15"
"set cmark 0"
"d b(t=49)"

"set string 1 c 15 0"
"set strsiz 0.3"
"draw string 5.65 7.9 "run.n""

"set string 1 c 15 90"
"set strsiz 0.2"
"draw string 0.6 4.4 height [km]"

"set string 1 c 15 0"
"set strsiz 0.2"
"draw string 5.65 0.7 buoyancy [m s`a-2`n]"

"printim ./figure/bpro_"run.n".png x2048 y1536"
"c"

"close 1"
n=n+1
endwhile

"quit"
