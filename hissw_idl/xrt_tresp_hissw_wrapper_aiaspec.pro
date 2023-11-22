



dates = {{ time }}
ifilters = {{ filters }}
xrt_path = {{ xrt_path }}
print, size(xrt_path)
print, 'Input Filters: ', ifilters

; Get the XRT temperature response for the below dates and save in format python likes
; This is an edited version to use the AIA default emission model (yes, AIA), recommended 
; by Kathy Reeves - see XRT Manual pg. 64 https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf
;
; 26-04-2021 IGH

; 12-10-2022 JMD EDIT to add Be-med and change date
; 5-9-2023 JMD CREATED THIS NEW VERSION (+edited to use AIA default emission model)
;~~~~~~~~~~~~~~~~~~~~~~~~~

print, dates
nd=n_elements(dates)

for d=0, nd-1 do begin
	aia=aia_get_response(/emiss, /full)
	modelname='AIA default spectrum'
	wave=aia.total.wave
	temp=10.^(aia.total.logte)
	spectrum = aia.total.emissivity
	abund_model= aia.general.abundfile
	ioneq_model = aia.general.ioneq_name
	dens_model = aia.general.model_name + ',p='+ trim(aia.general.model_pe,1)
	emiss_model = make_xrt_emiss_model(modelname, wave, temp, $
		spectrum, abund_model, ioneq_model, dens_model)

	wave_resp = make_xrt_wave_resp(contam_time=dates[d])
	;xrt_tresp = make_xrt_temp_resp(wave_resp, /apec_default)
	xrt_tresp = make_xrt_temp_resp(wave_resp, emiss_model)
	
	;print,wave_resp[1].contam.thick_time,wave_resp[1].contam.thick
	; For this only need 'Al-poly', 'Be-thin'
	filts=xrt_tresp.name
	for f=0, n_elements(filts)-1 do begin
	  idf=strpos(filts[f],';')
	  filts[f]=strmid(filts[f],0,idf)
	endfor
	
	ids=[]
	for ff=0, n_elements(ifilters)-1 do begin
		id=where(filts eq ifilters[ff])
		ids = [ids, id]
	endfor
	filters=filts[ids]
	
	print, 'Filters: ', filters

	print,xrt_tresp[ids].name
	units=xrt_tresp[ids[0]].temp_resp_units
	gd=where(xrt_tresp[ids[0]].temp gt 0.0)
	logt=alog10(xrt_tresp[ids[0]].temp[gd])
	tr=xrt_tresp[ids].temp_resp[gd]

	date=anytim(dates[d],/ccsds)
	
	t_in=dates[d]
	if n_elements(t_in) eq 0 then t_in = !stime
	extime = anytim(t_in,/ext)
	year = strtrim(string(extime(6)),1)
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
	

	save,file=xrt_path+'xrt_tresp_'+strmid(filname,0,8)+'_aiaspec.dat',filters,logt,tr,units,date
endfor

