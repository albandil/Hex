if ((GPVAL_VERSION == 4.3 || GPVAL_VERSION == 4.2) \
&&  (!strstrt(GPVAL_COMPILE_OPTIONS,"+OBJECTS"))) \
    print ">>> Skipping demo <<<\n" ; \
    print "This copy of gnuplot was built without support for placing rectangles\n" ; \
    exit ;

W = 0.5

set terminal png size 1024,400
set output "avail.png"

set multiplot

# compute number of lines in "avail.out"
num_lines = `cat avail.out | wc -l`

# set colour scheme
set palette model HSV functions 0.7,gray,1
set cbrange [0:7]

# set colours
do for [ll = 1:9] {
	set style line ll palette cb ll
}

# set frame
set logscale x
set xrange [0.1:100]
set yrange [0:2]
set grid
set xlabel "Ei - the impact energy in Rydbergs"
set ytics 10,10
set ylabel "initial state"
set cblabel "partial waves available"

spectr(li) = ((li==0) ? "s" : ((li==1) ? "p" : ((li==2) ? "d" : "?")))

set style rectangle fs solid noborder

do for [line_index = 1:num_lines] {

	# read line
	CMD_curr_line = sprintf("sed '%dq;d' avail.out", line_index)
	curr_line = system(CMD_curr_line)
	
	# get numbers
	ni = system(sprintf("echo \"%s\" | cut -f 1", curr_line))
	li = system(sprintf("echo \"%s\" | cut -f 2", curr_line))
	mi = system(sprintf("echo \"%s\" | cut -f 3", curr_line))
	minE = system(sprintf("echo \"%s\" | cut -f 4", curr_line))
	maxE = system(sprintf("echo \"%s\" | cut -f 5", curr_line))
	maxL = system(sprintf("echo \"%s\" | cut -f 6", curr_line))

	print "ni = ", ni
	print "li = ", li
	print "mi = ", mi
	print "minE = ", minE
	print "maxE = ", maxE
	print "maxL = ", maxL
	
	H = (ni-1)*(ni-1)+2*li+1+mi

	set ytics add ( sprintf("%s%s (m=%s)",ni,spectr(li),mi) H )

	# plot
#	plot (x < minE || x > maxE) ? 0/0 : H ls maxL+1 notitle
	set object line_index rect from minE,H-W to maxE,H+W fillcolor palette cb maxL
# 	plot line_index/10. lc palette cb line_index notitle
}

plot 0/0 palette cb 0 notitle

unset multiplot

set output
set terminal pop
