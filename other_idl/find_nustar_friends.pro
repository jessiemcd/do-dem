pro find_nustar_friends, doxrt=doxrt, doeis=doeis, downloadxrt=downloadxrt

;Finds data from other instruments that is coincident with NuTSTAR's pointing/timing during
;solar observations, and saves a bunch of handy information.
;
;Currently implemented for:
;
; -- Hinode/XRT - all filters, ignoring Ti-poly and C-poly after 14 June 2015
;
;=========================================================================================
;
;NuSTAR input: csv file with obsIDs, time intervals, and pointings. 

nuintervals = READ_CSV('solar_pointing_solarxy.csv', HEADER=nuhdr)

obsids = nuintervals.field02 ;nustar obsIDs for labeling, etc
tstart=convert_nustar_time(nuintervals.field05, /ut)
tstop=convert_nustar_time(nuintervals.field06, /ut)
nux = nuintervals.field09 ;should be in arcsec offset from solar center
nuy = nuintervals.field10 ;should be in arcsec offset from solar center


;=========================================================================================


coobs_stats=make_array(3,n_elements(obsids), value=0.)
coobs_stats[0,*]=obsids

filters=[]
for i=0, n_elements(obsids)-1 do begin
;for i=105, n_elements(obsids)-106 do begin
	keepdir = 0
	;Is there already a directory named for this obsID?
	obs = FILE_TEST(strtrim(obsids[i], 1)+'_coobs', /DIRECTORY)

	;print, obs

	;If not, let's make one!
	if obs EQ 0 then spawn, 'mkdir '+strtrim(obsids[i], 1)+'_coobs'

	;If there is one already, let's make sure we don't delete it later, even if we don't
	;end up adding anything new. 
	if obs EQ 1 then keepdir=1 
	
	;If the directory already exists, DON'T re-do!
	if obs EQ 1 then begin
		print, 'already did ', strtrim(obsids[i], 1)+'_coobs'
		CONTINUE
	endif
	
	;=============================================================
	;FIRST, LET'S CHECK FOR XRT COVERAGE
	;=============================================================
	print, [tstart[i], tstop[i]]
	print, [nux[i],nuy[i]]

	isthereany=0
	
	if doxrt EQ 1 then begin
	
		res = lets_get_that_xrt(nutimes=[tstart[i], tstop[i]], nupointing=[nux[i],nuy[i]], obsid=obsids[i], $
					checkimage=1)
					
				
		if res EQ !NULL and keepdir EQ 0 then spawn, 'rmdir '+strtrim(obsids[i], 1)+'_coobs'

		juf='-'
		if res NE !NULL then begin
			coobs_stats[1,i]=n_elements(res)
			f1=res.EC_FW1_
			f2=res.EC_FW2_
			filterlist=[]
			for j=0, n_elements(res)-1 do begin
				filterstring = f1[j]+'-'+f2[j]
				filterlist = [filterlist, filterstring]
			endfor
			uf = unique(filterlist, /sort)
			juf = uf.Join(', ')
			
			print, filterlist
			print, juf
			isthereany=1
		endif

		filters=[[filters], [juf]]
		
	endif
	
	if downloadxrt EQ 1 then begin

		if isthereany EQ 1 then begin

			res = lets_get_that_xrt(nutimes=[tstart[i], tstop[i]], nupointing=[nux[i],nuy[i]], obsid=obsids[i], $
					download=1, prep=1)
		endif

	endif
		
		
	
	;=============================================================
	;NOW LETS CHECK FOR EIS COVERAGE
	;=============================================================
	
	if doeis EQ 1 then begin
	
		eisres = lets_get_that_eis(nutimes=[tstart[i], tstop[i]], nupointing=[nux[i],nuy[i]], obsid=obsids[i], $
					checkimage=1)
		
		if eisres EQ !NULL and keepdir EQ 0 then spawn, 'rmdir '+strtrim(obsids[i], 1)+'_coobs'

				
		juf='-'
		if eisres NE !NULL then begin
			coobs_stats[2,i]=n_elements(eisres)
			
			;WHAT INFORMATION DO WE WANT TO KEEP ABOUT EIS?

		endif


	endif
	
	;=============================================================
	

endfor

if doxrt EQ 1 then begin
	write_csv, 'NuSTAR_XRT_coobs.csv', obsids, tstart, $
							tstop, nux, nuy, reform(fix(coobs_stats[1,*])), reform(filters), $
							Header=['obsID', 'Start', 'Stop', 'X Pointing', 'Y Pointing', '# Images', 'Which Filters']
endif


if doeis EQ 1 then begin
	write_csv, 'NuSTAR_EIS_coobs.csv', obsids, tstart, $
							tstop, nux, nuy, reform(fix(coobs_stats[2,*])), $
							Header=['obsID', 'Start', 'Stop', 'X Pointing', 'Y Pointing', '# EIS Images']
endif

stop
end