;RK 20.10.2017
;plot the baseline and line derived by
;fit_atran_base_and_line_parallel.pro
;RK 10.8.2018
;adding optional transmission (unsmoothed)
pro plot_bal_fit,prefix_fit,x,y,xrange=xrange,restframe=restframe,$
                 _extra=extra,overplot=overplot,points=points,resolution=resolution,$
                 transmission=transmission,tr_color=tr_color,params=params,$
                 in_wavelengths=in_wavelengths,wavelengths=wavelengths

  if ~(keyword_set(overplot) or $ 
       keyword_set(xrange) or $
       (n_elements(in_wavelengths) ge 2) )  then begin
     message,"Set at least overplot, xrange, or in_wavelengths",/info
     return
  endif
  if ~keyword_set(points) then points=500
  
  fns_fit=file_search('.',prefix_fit+'_*.fits',count=count)
  if count ne 10 then begin
     message,/info,"Fit files not found"
     return
  endif


  
  params=(readfits(fns_fit[0],hd))[x,y]
  get_za_alt,hd,za,alt
  
  for i=1,9 do begin
     if (i ne 1) and (i ne 7) then begin
        param=(readfits(fns_fit[i]))[x,y]
        params=[[[[params]]],[[param]]]
     endif
  endfor
  
  i_alt=0
;  i_bal=1
  i_center=2-1
  i_cont=3-1
  i_flux=4-1
  i_pwv=5-1
  i_resolution=6-1
;  i_residuals=7
  i_width=8-2
  i_zero=9-2
  
  if ~keyword_set(restframe) then begin
     z=sxpar(hd,'REDSHIFT')
     params[i_center]*=z+1.
     params[i_width]*=z+1.
  endif

;resorting paramas in the right order
  if keyword_set(resolution) then begin
     params=[params[[i_alt,i_pwv]],za,params[[i_cont,i_zero,i_flux,i_center,i_width]],resolution]
  endif else begin
     params=[params[[i_alt,i_pwv]],za,params[[i_cont,i_zero,i_flux,i_center,i_width,i_resolution]]]
  endelse

if n_elements(in_wavelengths) lt 2 then begin  
   if keyword_set(overplot) then begin
      if !X.type eq 0 then begin
         wavelengths=!X.crange[0]+(!X.crange[1]-!X.crange[0])*findgen(points)/(points-1.)
      endif else begin
         wavelengths=10^(!X.crange[0]+(!X.crange[1]-!X.crange[0])*findgen(points)/(points-1.))
      endelse
   endif else begin     
      if ~keyword_Set(xtype) then begin
         wavelengths=Xrange[0]+(Xrange[1]-Xrange[0])*findgen(points)/(points-1.)
      endif else begin
         wavelengths=10^(Xrange[0]+(Xrange[1]-Xrange[0])*findgen(points)/(points-1.))
      endelse
   endelse
endif else begin
   wavelengths=in_wavelengths
endelse

  if keyword_set(overplot) then begin
     oplot,wavelengths,atran_base_and_line(wavelengths,params),_extra=extra
  endif else begin
     plot,wavelengths,atran_base_and_line(wavelengths,params),_extra=extra
  endelse
  oplot,wavelengths,atran_base_and_line(wavelengths,params*[1,1,1,1,1,0.,1,1,1]),linest=1
  
  if keyword_set(transmission) then begin
     if n_elements(tr_color) eq 0 then tr_color=cgcolor('hot pink')
     axis,/yaxis,yr=[0,1],/save,color=tr_color,ytitle='atmospheric transmission'
     oplot,wavelengths,atran_transmission(wavelengths,alt,params[1],za),color=tr_color
  endif
end
