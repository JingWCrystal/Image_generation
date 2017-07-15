WIDTH=100
HEIGHT=100

scan_left=int(bp-WIDTH/2)
scan_right=scan_left+WIDTH

read_package=[]
for read in samfile.fetch(scan_left,scan_right):
	if read_can_shown(read,scan_left,scan_right):
		new_bases,new_base_quality=rearrange_string(read)
		read_package.append((new_bases,new_base_quality,read))
dic=read_to_dictionary(read_package, scan_right,HEIGHT)
draw_pgn('left',dic,width,height,scan_left,scan_right,img_name)

# for right breakpoint, 
# draw_pgn('right',dic,width,height,scan_left,scan_right,img_name)

def read_can_shown(read,scan_left,scan_right):
	if (read.cigarstring != None) and (read.start >=scan_left) and (read.end<=scan_right) and (read.mapping_quality >10) :
		return ture
	else:
		return false

def rearrange_string(read):
	bases=read.sequence   
	new_bases=''
	new_base_quality=[]

	read_pos=0
	for a_cigar in read:
		if a_cigar=='M':
			new_bases=new_bases+ bases[read_pos]
			new_base_quality.append( read.query_qualities[read_pos])
			read_pos+=1
		elif a_cigar=='I': 
			read_pos+=1
		elif a_cigar=='D': # the deletion partion is not calculated in "read length"
			new_bases=new_bases+'d'	
			new_base_quality.append(0)
		elif a_cigar=='S':
			new_bases=new_bases+'s'	
			new_base_quality.append(-1)
			read_pos+=1
	
    return new_bases,new_base_quality       

def read_to_dictionary(read_package, scan_right,HEIGHT):
	dictionary={}
	read_list=[0 for i in range(HEIGHT)] 
	for i in range(HEIGHT):
		read_list[i]=[]
	row_ptr=0 
	for base_and_read in read_package:
		if row_ptr<HEIGHT:
			base=base_and_read[0]
			quality=base_and_read[1]
			read=base_and_read[2]
			if read.cigartuples[0][0]==4 : 
				read_pos1=read.start - read.cigartuples[0][1]
				read_pos2=read_pos1+read_length
			else:
				read_pos1=read.start
				read_pos2=read_pos1+read_length
			if read.is_paired:  # read not paired，do not consider this kind of read
				is_concordant=read.is_proper_pair
				if 'S' in read.cigarstring:
					is_clipped=True
				else:
					is_clipped=False
				if is_empty(read_list):
					read_list[row_ptr].append( (read_pos1,read_pos2,base,quality,is_clipped,is_concordant) )
				else: 
					row_ptr = get_shortest_tail_row(read_list,scan_r_pos)
					if read_pos1 >= read_list[row_ptr][-1][1]: 
						read_list[row_ptr].append( (read_pos1,read_pos2,base,quality,is_clipped,is_concordant) )
					else:  
						row_ptr=find_next_empty_row(read_list)
						read_list[row_ptr].append( (read_pos1,read_pos2,base,quality,is_clipped,is_concordant))

	for i in range(HEIGHT):
		dictionary[i]=read_list[i]
	return dictionary

def draw_pgn(which_bp,dic,width,height,scan_left,scan_right,img_name):
	newIm = Image.new ("RGB", (width,height),(255,255,255))
	for key in range(height):
		for read_tuple in dic[key]:
			read_pos1=read_tuple[0]
			read_pos2=read_tuple[1]
			base=read_tuple[2]
			quality=read_tuple[3]
			is_clipped=read_tuple[4]
			is_concordant=read_tuple[5]
			col=read_pos1-scan_left
			index_in_read=0
			for i in range(len(base)):
				if col >=0  and col < width :
					row=key
					# print row,col
					red,green,blue=get_RGB(which_bp,base[index_in_read],quality[index_in_read],is_clipped,is_concordant)
					newIm.putpixel((col,row),(red,green,blue))
					index_in_read=index_in_read+1
					col=col+1
				elif col<0:
					index_in_read=index_in_read+1
					col=col+1
	newIm.save(img_name,"PNG")

def get_RGB(which_bp,base,quality,is_clipped,is_concordant):
	if is_clipped :
		red=255
		green=0
		blue=0
		if quality>1:
			green=green+255-6*quality
			blue=blue+255-6*quality
			return red, green, blue
		elif quality== 0 or quality== -1 :
			return 255,255,255
	elif is_concordant:
		red=0
		green=255
		blue=0
		if quality>1:
			red=red+255-6*quality
			blue=blue+255-6*quality
			return red, green, blue
		elif quality==0 or quality== -1: 
			return 255,255,255
	elif not is_concordant:
		red=0
		green=0
		blue=255
		if quality>1:
			red=red+255-6*quality
			green=green+255-6*quality
			return red, green, blue
		elif quality==0 or quality== -1:
			return 255,255,255