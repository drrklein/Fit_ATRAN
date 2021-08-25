;RK 23.3.2017
;load cube and (o)plot the spectrum for spaxel (x,y)
;continuum can only be reconstructed for constant and linear baselines

pro plot_spaxel,fn_cube,x,y,restframe=restframe,err_ext=err_ext,$
                cube_ext=cube_ext,wave_ext=wave_ext,wavehdr=wavehdr,_extra=extra,overplot=overplot
  if n_elements(cube_ext) eq 0 then cube_ext=3
  if n_elements(wave_ext) eq 0 then wave_ext=5
  
  
  cube=mrdfits(fn_cube,cube_ext)

  if ~keyword_set(wavehdr) then begin
     wave=mrdfits(fn_cube,wave_ext)
  endif else begin
     wave=wave_hdr(headfits(fn_cube))
  endelse
  

  if keyword_set(restframe) then begin
     dummy=readfits(fn_cube,hd)
     z=sxpar(hd,'REDSHIFT')
     wave/=1.+z
  endif
  
  if keyword_set(overplot) then begin
     oplot,wave,cube[x,y,*],_extra=extra
  endif else begin
     plot,wave,cube[x,y,*],_extra=extra
  endelse

  if n_elements(err_ext) eq 1 then begin
     err=mrdfits(fn_cube,err_ext)
     oplotyerr,wave,cube[x,y,*],err[x,y,*],_extra=extra
  endif
  
end
