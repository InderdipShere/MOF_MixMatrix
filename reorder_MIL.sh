# change the xaxis of the plot
file=$1
box=$2
cen=$3
awk '{if($1=="#")
	    print "#"
	  else if($1==None)
	      print "  " 
	    else {
	      {abox=(($1-'$cen')/'$box')}
	      if (abox>=0.5)
            {abox=1}
            print  ($1-'$cen')-(int(abox)*'$box'), $2 }}' $file
#print ($1-'$cen')-(int(($1-'$cen')/'$box')*'$box'), $2}' $file
