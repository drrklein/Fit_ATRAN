;RK 9.10.2016
; the ATRAN files have different wavelength grids. Regridding allows
; to save a lot of memorry

if n_elements(atran) ne 1800 then restore,'atran_fifils.sav',/ver

alt=38+lindgen(8)
pwv=[1,2,3,4,5,6,7,8,10,12,14,15,16,18,20,22,25,27,30,32,35,37,40,45,50]
ZA=30+5*lindgen(9)


dummy_hd=["CHANNEL = 'BLUE    '           / Detector channel",$
          "G_ORD_B =                    2 / Blue grating order to be used"]
order = 2
ch='b'

min=50.
max=205.

over_sample=10^0.75
wave_grid=[min]
repeat begin
   R=fifi_ls_get_resolution(dummy_hd,wmean=wave_grid[-1])
   if (order eq 2) and (wave_grid[-1] ge 70) then begin
      dummy_hd[1]="G_ORD_B =                    1 / Blue grating order to be used"
      order = 1
   endif
   if (ch eq 'b') and (wave_grid[-1] ge 118) then begin
      dummy_hd[0]="CHANNEL = 'RED     '           / Detector channel"
      ch='r'
   endif
   wave_grid=[wave_grid,wave_grid[-1]*(1.+1./(over_sample*R))]
endrep until wave_grid[-1] ge max

for a=0,n_elements(alt)-1 do for z=0,n_elements(ZA)-1 do for p=0,n_elements(pwv)-1 do begin
   fn=string(alt[a],ZA[z],pwv[p],format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
   trans=interpol(atran[fn].trans,atran[fn].wave,wave_grid)
   trans=gauss_smooth(trans,over_sample/2.3458); R seems to be lambda/FWHM ; here sigma is needed
   atran[fn]=trans
endfor

save,atran,wave_grid,fn,alt,pwv,za,filename='atran_fifils_regrid.sav'



end
