ADD_PATH, '~/ssw/gen/idl/util/'
ADD_PATH, '~/ssw/sdo/aia/idl/response/'


;  tresp=aia_get_response(/temperature,/dn,/eve,timedepend_date='01-Jul-2010')
tresp=aia_get_response(/temperature,/dn,/eve) 
ids=[0,1,2,3,4,6]
channels=tresp.channels[ids]
logt=tresp.logte
tr=tresp.all[*,ids]
units=tresp.units
save,file='aia_tresp_en.dat',channels,logt,tr,units 

; Do the noblend version so don't include crosstalk between 131+335 and 94+304
tresp=aia_get_response(/temperature,/dn,/eve,/noblend)
ids=[0,1,2,3,4,6]
channels=tresp.channels[ids]
logt=tresp.logte
tr=tresp.all[*,ids]
units=tresp.units
save,file='aia_tresp_en_nb.dat',channels,logt,tr,units
