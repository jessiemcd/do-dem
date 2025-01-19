FUNCTION lets_get_that_eis, nutimes=nutimes, nupointing=nupointing, obsid=obsid, checkimage=checkimage

	; Written by Jessie Duncan, September 2022
	;-Developing an automated search process for getting NuSTAR/EIS co-observations
	;-Assuming NuSTAR search is handled elsewhere, so here we will start with NuSTAR times and 
	;pointings.

	;Keywords: 
	;			nutimes –– nustar data interval (or other interesting interval of time), see default for
	;						format. 
	;			nupointing –– nustar pointing (or other coordinate of interest)
	;			obsid –– nustar obsid. only used to print on the resulting pointing image, so optional
	

	;default times and pointing are from Hi-C co-observation nustar campaign 
	default, nutimes, '2018-may-29 '+['1500','2100']
	default, nupointing, [-110., 250]
	default, checkimage, 0

	default, obsid, 000
	
	spawn, 'mkdir '+strtrim(obsid, 1)+'_coobs'
	thepath = './'+strtrim(obsid, 1)+'_coobs/'
	print, thepath
	
	nuFOV = [10., 10.] ;Approximate "guaranteed" NuSTAR FOV in arcminutes, assume centered around pointing 
	nuextent = [5.,5.,5.,5.]*60. ;NuSTAR "guaranteed" FOV EXTENT (in arcseconds) in each direction from the pointing
	nuextent = [6.,6.,6.,6.]*60. ;NuSTAR "guaranteed" FOV EXTENT (in arcseconds) in each direction from the pointing
	xmin = nupointing[0]-nuextent[0]
	xmax = nupointing[0]+nuextent[1]
	ymin = nupointing[1]-nuextent[2]
	ymax = nupointing[1]+nuextent[3]
	
	s=eis_obs_structure(nutimes[0],nutimes[1])
	
	if typename(s) EQ 'INT' then begin
		print, 's int'
		return, !NULL
	endif
	
	
	;THIS WOULD NOT CATCH A SINGLE FILE, FIX IN FUTURE
	if typename(s) EQ 'ANONYMOUS' or typename(s) EQ 'STRUCT' then s2=eis_filter_obs_struc(s, xcen=[xmin,xmax],ycen=[ymin,ymax])
	
	
	if typename(s2) EQ 'INT' then begin
		print, 's2 int'
		return, !NULL
	endif
	
	
	if checkimage EQ 1 then begin

			;Makes plot of Sun with full-disk image, and FOV boxes for each EIS file in the time interval: 
			;Green boxes for FOVs containing NuSTAR pointing
			;Red boxes for FOVs NOT containing NuSTAR pointing
			;Pink X shows input NuSTAR pointing on solar disk
	
	

			;Find a full disk AIA image

			searchfile = vso_search(nutimes[0],nutimes[1], instr='aia', wave='94')
			getfile = vso_get(searchfile[0],filenames=fulldiskcontext)
			aia_lct, r,g,b,wave=94, /load

	
			spawn, 'mv '+file_basename(fulldiskcontext)+' ./'+thepath
			fits2map, thepath+file_basename(fulldiskcontext), exmap
			
			if exmap EQ !NULL then begin
				print, ""
				print, 'AIA Download Issue!'
				print, nutimes[0]
				print, ""
				print, ""
				print, ""
				return, s2
			endif
			
			
			popen, thepath+'eis_coobs_'+strtrim(obsid,1)+'.ps', xsi=7, ysi=10
			!Y.margin=4.
			!X.margin=4.
			ch=1.1
			th=6
			lnth=4
			!p.multi=[0,1,1]

			th=7
			fth=4
			;===============================================================================
			;===============================================================================

			
			reverse_ct
			plot_map, exmap, /limb, xrange=[-1500, 1500], yrange=[-1500, 1500], /log, col=255, lcol=255
			hsi_linecolors
			oplot, [nupointing[0], nupointing[0]], [nupointing[1], nupointing[1]], psym=7, thick=3, symsize=2, col=2
			
			;Make BOXES!
			if s NE !NULL then begin
				for i=0, n_elements(s)-1 do begin
					;print, i
					fovx = s[i].FOVX
					fovy = s[i].FOVY					
					target = [s[i].XCEN, s[i].YCEN]

					box = fltarr(2,5)
					box[0,*] = fovx*0.5*[-1,-1,1,1,-1]
					box[1,*] = fovy*0.5*[1,-1,-1,1,1]


					theta=0
					cosx = cos(theta*!pi/180.)
					sinx = sin(theta*!pi/180.)
					rotarr = [[cosx,-sinx],[sinx,cosx]]
					for j=0,4 do box[*,j] = box[*,j]#rotarr + target
	
					oplot, box[0,*], box[1,*], col=6, thick=5
				endfor
			endif
	
			;Make BOXES!
			if s2 NE !NULL then begin
				for i=0, n_elements(s2)-1 do begin
					;print, i
					fovx = s2[i].FOVX
					fovy = s2[i].FOVY
					target = [s2[i].XCEN, s2[i].YCEN]

					box = fltarr(2,5)
					box[0,*] = fovx*0.5*[-1,-1,1,1,-1]
					box[1,*] = fovy*0.5*[1,-1,-1,1,1]


					theta=0
					cosx = cos(theta*!pi/180.)
					sinx = sin(theta*!pi/180.)
					rotarr = [[cosx,-sinx],[sinx,cosx]]
					for j=0,4 do box[*,j] = box[*,j]#rotarr + target
	
					oplot, box[0,*], box[1,*], col=3, thick=5
				endfor
			endif
			
		charsize=1
			charth=1.5
			xyouts, -1200.,-1200., 'Approx. NuSTAR pointing from '+strtrim(obsid,1), charsize=charsize, charth=charth
			xyouts, -1200.,-1300., 'Interval: '+strtrim(nutimes[0], 1)+' to '+strtrim(nutimes[1], 1), charsize=charsize, charth=charth
			xyouts, -1200., -1400., strtrim(n_elements(s2), 1)+' total EIS files in interval', charsize=charsize, charth=charth
	
			;===============================================================================
			;===============================================================================

			!p.multi=0
			!Y.margin=[4.,2.]
			pclose
			;spawn, 'open '+thepath+'eis_coobs_'+strtrim(obsid,1)+'.ps'


			;===============================================================================
			;===============================================================================


	endif


print, s2

return, s2



;stop
end