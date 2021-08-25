;RK 29.3.2015
;create the wavelengths from FITS header keywords

function wave_hdr,hd,naxis3=naxis3,crpix3=crpix3,cdelt3=cdelt3,crval3=crval3,$
                  _extra=extra
  naxis3=sxpar(hd,'NAXIS3',_extra=extra)
  crpix3=sxpar(hd,'CRPIX3',_extra=extra)
  cdelt3=sxpar(hd,'CDELT3',_extra=extra)
  crval3=sxpar(hd,'CRVAL3',_extra=extra)
  wave=(findgen(naxis3)-(crpix3-1))*cdelt3+crval3
 
  return,wave
end
