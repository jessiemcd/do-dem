;+
; Name: hsi_linecolors
;
; Purpose:  Define discrete colors to use for HESSI plots
;
; Calling Sequence:  hsi_linecolors
;
; Input Keywords:
;  pastels - if set, add some pastel colors
;  help - if set, show a colorbar showing indices of new colors
;  
; Output: The color table will be changed
; 
; Examples:
;   hsi_linecolors
;   hsi_linecolors, /help, /pastel
;   hsi_linecolors,/pastel
;   
; Method:  Get current color table, redefine first 11 indices
;	of r,g,b arrays and reload color table.
;   Result is that setting color to any of the first 11 indices
;	will give the following color:
;	0 - black
;	1 - white
;	2 - magenta
;	3 - green
;	4 - cyan
;	5 - yellow
;	6 - red
;	7 - blue
;	8 - orange
;	9 - olive
;	10 - purple
;	11 - violet
;	12 - dark green
;	13 - light red
;	14 - light green
;	15 - light blue
;	16 - light yello
;	17 - light purple
;	18 - light gray
;	19 - gray
;
;Written:  Kim Tolbert, 25- Mar-2002
;Modifications:  
; 11-Nov-2003, Kim.  Added violet and dark green
; 31-May-2012, Kim.  Added pastel colors via /pastel option
;
;-

pro hsi_linecolors, pastels=pastels, help=help

tvlct, r,g,b, /get

r[0:12] = [0, 255, 255,   0,   0, 213, 255,   0, 255, 128, 128, 128,   0]
g[0:12] = [0, 255,   0, 255, 255, 213,   0,   0, 128, 128,   0,   0, 128]
b[0:12] = [0, 255, 255,   0, 255,   0,   0, 255,   0,   0, 128, 170,   0]
used = 13

if keyword_set(pastels) then begin
;  r[13:19] = [255, 240, 230, 255, 240, 240, 150]
;  g[13:19] = [240, 255, 250, 255, 240, 240, 150]
;  b[13:19] = [240, 240, 255, 240, 255, 240, 150]
  r[13:19] = [255, 229, 230, 255, 230, 240, 150]
  g[13:19] = [240, 244, 250, 255, 230, 240, 150]
  b[13:19] = [240, 229, 255, 240, 245, 240, 150]
  used = 20
endif

tvlct, r,g,b

if keyword_set(help) and ( (!d.name eq 'X') or (!d.name eq 'WIN') ) then begin
   bar=rebin(indgen(used),32*used,used*4,/sample)
   wtemp=!d.window
   wdef,zz,/ur,image=bar
   tv,bar
   xyouts,indgen(used)*32+8, 32,strtrim(sindgen(used),2),/device,size=1.5,charthick=2.
   if wtemp ne -1 then wset,wtemp
endif

end
