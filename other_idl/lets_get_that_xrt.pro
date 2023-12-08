
function break_time, t_in, qdebug = qdebug
  if n_elements(t_in) eq 0 then t_in = !stime
  extime = anytim(t_in,/ext)
  year = strtrim(string(extime(6)),1)
  month = strtrim(string(extime(5)),1)
  if strlen(month) eq 1 then month = '0'+month
  day = strtrim(string(extime(4)),1)
  if strlen(day) eq 1 then day = '0'+day
  hour = strtrim(string(extime(0)),1)
  if strlen(hour) eq 1 then hour = '0'+hour
  minute = strtrim(string(extime(1)),1)
  if strlen(minute) eq 1 then minute = '0'+minute
  second = strtrim(string(extime(2)),1)
  if strlen(second) eq 1 then second = '0'+second
  filname = year+month+day+'_'+hour+minute+second
  if keyword_set(qdebug) then stop
  return, filname
end





FUNCTION lets_get_that_xrt, nutimes=nutimes, nupointing=nupointing, obsid=obsid, $
				download=download, checkimage=checkimage, prep=prep, path=path

		; Written by Jessie Duncan, September 2022
		;-Developing an automated search process for getting NuSTAR/XRT co-observations
		;-Assuming NuSTAR search is handled elsewhere, so here we will start with NuSTAR times and 
		;pointings.
		
		
		; Note: you will need some of the database files available via SSWDB: 
		;https://www.mssl.ucl.ac.uk/surf/sswdoc/solarsoft/sswdb_install.html
		;Particularly, xrt_genxcat_sirius + xrt_genxcat
		;If your installation isn't working, check to make sure you're using the 'spawn' keyword
		; i.e. sswdb_upgrade, /spawn 
		; after downloading the db file and placing it in $SSW/site/setup/

		;Keywords: 
		;			nutimes –– nustar data interval (or other interesting interval of time), see default for
		;						format. 
		;			nupointing –– nustar pointing (or other coordinate of interest) in x,y arcseconds from 
		;                       solar center
		;			obsid –– nustar obsid. only used to print on the resulting pointing image, so optional
		;			download -- set download=1 to download all the XRT files in the time range (note: currently
		;						before any pointing cuts)
		;           prep -- If downloading, also run xrt_prep on the data files.
		;			checkimage -- set checkimage=1 to plot a full-disk solar image, with the NuSTAR pointing as 
		;						a green X and the FOV of every XRT image shown as red boxes
		;
		
		
;                Sample run: 
;                
;                IDL> nutimes='2018-may-29 '+['1909','1955'] 
;                IDL> nupointing=[-110., 250]
;                IDL> lgtx = lets_get_that_xrt(nutimes=nutimes, nupointing=nupointing, $
;                IDL>    download=1, prep=1, obsid=80410203001)
;                
;                This will download level 0 files, and place them in './80410203001_coobs/'
;                It will then prep files, and place level 1 files + grade maps in './80410203001_coobs/XRT_for_DEM/'	
;			

		default, download, 0
		default, checkimage, 0
		;default times and pointing are from Hi-C co-observation nustar campaign 
		default, nutimes, '2018-may-29 '+['1900','2000']
		default, nupointing, [-110., 250]

		default, obsid, 000
		default, prep, 0


		thepath = './'+strtrim(obsid, 1)+'_coobs/'
		print, thepath
		
		if keyword_set(path) then begin
			thepath = path
		endif
		
		res = FILE_TEST(thepath, /DIRECTORY)
		if res EQ 0 then begin
			FILE_MKDIR, thepath
		endif


		;DO NOT USE WITHOUT ALSO CHANGING THE OBSID AND THE POINTING!
		;ishikawa and krucker orbit 1
		;nutimes = '2010-jul-13 '+['0531','0549']
		;orbit 2
		;nutimes = '2010-jul-13 '+['0856','0930']
		;orbit 3
		;nutimes = '2010-jul-13 '+['1330','1419']
		;orbit 4
		;nutimes = '2010-jul-13 '+['1719','1741']

		;foxsi-2
		;nutimes = '2014-dec-11 '+['1908','1914'] 
		;Hi-C NuSTAR Co-observation
		;nutimes = '2018-may-29 '+['1700','1800']

		;nupointing = [-110., 250]
		;nupointing = [-450., 250]



		nuFOV = [10., 10.] ;Approximate "guaranteed" NuSTAR FOV in arcminutes, assume centered around pointing 
		nuextent = [5.,5.,5.,5.]*60. ;NuSTAR "guaranteed" FOV EXTENT (in arcseconds) in each direction from the pointing
		xmin = nupointing[0]-nuextent[0]
		xmax = nupointing[0]+nuextent[1]
		ymin = nupointing[1]-nuextent[2]
		ymax = nupointing[1]+nuextent[3]


		xrt_cat, nutimes[0], nutimes[1], catx, ofiles, urls=1

		if n_elements(catx) EQ 1 then begin
			print, 'Nothing Found in this Time Interval'
			return, []
			;stop
		endif


		if download EQ 1 then begin
			for i=0, n_elements(ofiles)-1 do begin
				filestring = string(ofiles[i])
				filename = file_basename(ofiles[i])
				w = wget(filestring, filename=filename, directory=thepath)
				print, w
			endfor
			
			if prep EQ 1 then begin
				xrt_files = file_search(thepath+'XRT*.fits')
				print, xrt_files
	
				for i=0, n_elements(xrt_files)-1 do begin
					read_xrt,xrt_files[i],ind,data,/force
					xrt_prep,ind,data,indp,datap,/float,grade_map=gm,grade_type=0,/coalign
					indp.timesys='UTC'
					filtout=indp.ec_fw1_+'_'+indp.ec_fw2_
					resout=strcompress(string(indp.naxis1),/rem)
					fnameout='XRT_'+break_time(indp.date_obs)+'_'+filtout+'_'+resout+'.fits'
					write_xrt,indp,datap,outdir=thepath+'/XRT_for_DEM',outfile=fnameout,/ver, /MAKE_DIR
					gfnameout='gm_XRT_'+break_time(indp.date_obs)+'_'+filtout+'_'+resout+'.fits'
					write_xrt,indp,gm,outdir=thepath+'/XRT_for_DEM',outfile=gfnameout,/ver, /MAKE_DIR
				endfor
			endif
		endif


		;Deprecated option: Finding files where 
		;the western edge of the XRT FOV is west of the NuSTAR x pointing 
		; 						AND
		;the eastern edge of the XRT FOV is east of the NuSTAR x pointing
		;print, 'Western XRT X boundaries:'
		;print, catx.XCEN+catx.FOVX/2
		;print, 'Eastern XRT X boundaries:'
		;print, catx.XCEN-catx.FOVX/2
		;print, 'XRT naxis:'
		;print, catx.naxis1
		;print, ''

		print, 'Initial List:'
		print, size(catx)

		;Finding files where the x coordinate of the nustar pointing is between the eastern and western edges of the XRT FOV.
		xfiles = catx[where(catx.XCEN+catx.FOVX/2. GE nupointing[0] AND catx.XCEN-catx.FOVX/2. LE nupointing[0], complement=nonx, /null)]
		;Within the above, finding files where the y coordinate of the nustar pointing is between the northern and southern edges
		;of the XRT FOV.
		files=[]
		if xfiles NE !null then files = xfiles[where(xfiles.YCEN+xfiles.FOVY/2. GE nupointing[1] AND xfiles.YCEN-xfiles.FOVY/2. LE nupointing[1], /null)]

		;Doing the same selection but for the list of file URLs
		xlinks = ofiles[where(catx.XCEN+catx.FOVX/2. GE nupointing[0] AND catx.XCEN-catx.FOVX/2. LE nupointing[0], /null)]
		links=[]
		if xlinks NE !null then links = xlinks[where(xfiles.YCEN+xfiles.FOVY/2. GE nupointing[1] AND xfiles.YCEN-xfiles.FOVY/2. LE nupointing[1], /null)]

		;Defining files where the X coordinate of the nustar pointing is outside of the XRT FOV x range, using the complement of the 
		;earlier where statement
		;Defining non-files as null to start, will be updated IFF there are any.
		nonfiles=[]
		if nonx NE !null then begin
			nonxfiles = catx[nonx]
			;Defining where the "non-x" files are IN the XRT FOV y range (complement will be when they are NOT)
			tt = nonxfiles[where(nonxfiles.YCEN+nonxfiles.FOVY/2. GE nupointing[1] AND nonxfiles.YCEN-nonxfiles.FOVY/2. LE nupointing[1], complement=nony, /null)]
			;If there are any in the complement, define the nonfiles.
			if nony NE !null then nonfiles = nonxfiles[nony]
		endif
		
		if n_elements(files) EQ 0 then return, []
		
		;REMOVING TI-POLY AND C-POLY RESULTS (LIGHT LEAKS)
		files = files[where(files.EC_FW1_ NE 'C_poly' and files.EC_FW2_ NE 'Ti_poly', /null)]
		
		if n_elements(files) EQ 0 then return, []
		
		;REMOVING SUB-1s DATA CAPTURES
		files = files[where(files.EXPTIME GT 1., /null)]
		
		if n_elements(files) EQ 0 then return, []

		print, 'X List:'
		print, size(xfiles)
		print, 'Full List:'
		print, size(files)
		print, 'No X List:'
		print, size(nonxfiles)
		print, 'No List:'
		print, size(nonfiles)
		print, ''

		
		print, 'Files with FOV Match - filter combinations:'

		for i=0, n_elements(files)-1 do begin
			print, files[i].EC_FW1_, '  ', files[i].EC_FW2_, '  ', files[i].EXPTIME, '  ', $
					files[i].FOVX, '  ', files[i].FOVY
		endfor

		write_csv, thepath+strtrim(obsid,1)+'_XRT_coobservations.csv', files.EC_FW1_, files.EC_FW2_, $
						files.EXPTIME, files.FOVX, files.FOVY, $
						Header=['Filter 1', 'Filter 2', 'Exposure Time (s)', 'FOV Size X', 'FOV Size Y']
		coobs_XRT_headers = files
		SAVE,coobs_XRT_headers, FILENAME = thepath+strtrim(obsid,1)+'_XRT_coobs.sav'
		SAVE, links, FILENAME = thepath+strtrim(obsid,1)+'_XRT_coobs_filelinks.sav'

		if checkimage EQ 1 then begin

			;Makes plot of Sun with full-disk image, and FOV boxes for each XRT file in the time interval: 
			;Green boxes for FOVs containing NuSTAR pointing
			;Red boxes for FOVs NOT containing NuSTAR pointing
			;Pink X shows input NuSTAR pointing on solar disk
	
	
			;Find a full-disk XRT image if available
			;fullfiles = where(catx.naxis1 eq 2048 and catx.naxis2 eq 2048, /null)
			fullfiles = where(catx.FOVX eq 2106.5701 and catx.FOVY eq 2106.5701, /null)
			if fullfiles NE !NULL then begin
				index=fullfiles[0]
				fulldiskcontext = ofiles[index]
				w=WGET(fulldiskcontext, FILENAME=file_basename(fulldiskcontext))
				loadct, 1
			;If there isn't one, find a full disk AIA image
			endif else begin
				print, 'You have no full disk XRT in this interval, downloading an AIA94 image'
				searchfile = vso_search(nutimes[0],nutimes[1], instr='aia', wave='94')
				getfile = vso_get(searchfile[0],filenames=fulldiskcontext)
				aia_lct, r,g,b,wave=94, /load
			endelse
	
			spawn, 'mv '+file_basename(fulldiskcontext)+' ./'+thepath
			fits2map, thepath+file_basename(fulldiskcontext), exmap
			
			if exmap EQ !NULL then begin
				print, ""
				print, 'Download Issue!'
				print, nutimes[0]
				print, ""
				print, ""
				print, ""
				return, files
			endif
			
			
			popen, thepath+'coobs_'+strtrim(obsid,1)+'.ps', xsi=7, ysi=10
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
			if files NE !NULL then begin
				for i=0, n_elements(files)-1 do begin
					;print, i
					fov = files[i].FOVX
					target = [files[i].XCEN, files[i].YCEN]

					box = fltarr(2,5)
					box[0,*] = fov*0.5*[-1,-1,1,1,-1]
					box[1,*] = fov*0.5*[1,-1,-1,1,1]


					theta=0
					cosx = cos(theta*!pi/180.)
					sinx = sin(theta*!pi/180.)
					rotarr = [[cosx,-sinx],[sinx,cosx]]
					for j=0,4 do box[*,j] = box[*,j]#rotarr + target
	
					oplot, box[0,*], box[1,*], col=3, thick=5
				endfor
			endif
			
			if nonfiles NE !NULL then begin
				for i=0, n_elements(nonfiles)-1 do begin
					;print, i
					fov = nonfiles[i].FOVX
					target = [nonfiles[i].XCEN, nonfiles[i].YCEN]

					box = fltarr(2,5)
					box[0,*] = fov*0.5*[-1,-1,1,1,-1]
					box[1,*] = fov*0.5*[1,-1,-1,1,1]


					theta=0
					cosx = cos(theta*!pi/180.)
					sinx = sin(theta*!pi/180.)
					rotarr = [[cosx,-sinx],[sinx,cosx]]
					for j=0,4 do box[*,j] = box[*,j]#rotarr + target
	
					oplot, box[0,*], box[1,*], col=6, thick=5
				endfor
			endif
	
			charsize=1
			charth=1.5
			xyouts, -1200.,-1200., 'Approx. NuSTAR pointing from '+strtrim(obsid,1), charsize=charsize, charth=charth
			xyouts, -1200.,-1300., 'Interval: '+strtrim(nutimes[0], 1)+' to '+strtrim(nutimes[1], 1), charsize=charsize, charth=charth
			xyouts, -1200., -1400., strtrim(n_elements(files), 1)+' Total XRT Images in Interval', charsize=charsize, charth=charth
	
			;===============================================================================
			;===============================================================================

			!p.multi=0
			!Y.margin=[4.,2.]
			pclose
			;spawn, 'open '+thepath+'coobs_'+strtrim(obsid,1)+'.ps'


			;===============================================================================
			;===============================================================================

		endif
return, files

;stop
end