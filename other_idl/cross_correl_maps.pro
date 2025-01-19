FUNCTION	cross_correl_maps, maps

	corrections = fltarr( 2, n_elements(maps) )
	for i=0, n_elements(maps)-2 do corrections[*,i:i+1] += get_correl_offsets( maps[i:i+1].data )
	corrections[0,*] *= maps[0].dx
	corrections[1,*] *= maps[0].dy

	; Those are the corrections between consecutive maps.  Now take the cumulative sum to 
	; get the corrections from (for now) the first map.
	for i=n_elements(maps)-1, 0, -1 do corrections[0,i] = total(corrections[0,0:i])
	for i=n_elements(maps)-1, 0, -1 do corrections[1,i] = total(corrections[1,0:i])

	; Apply the alignment to the maps.
	selfcor = maps
	for i=0, n_elements(maps)-1 do selfcor[i] = shift_map( maps[i], corrections[0,i], corrections[1,i] )

	return, selfcor

END
