;RK 7.6.2019
; coding gauss_smooth so that it works for sigmas > 100

function mygauss_smooth,data,sigma_ext,kernel=kernel,width=width,_extra=extra
  sigma=sigma_ext
  n_dim=size(data,/n_dim)
  n_ele=n_elements(sigma)

  if n_dim gt 8 then message,'Must be less than 8 dimesions'
  if n_dim ne n_ele then begin
     if n_ele ne 1 then begin
        message,'Number of dimensions and elements of sigma must match or sigma must be a scalar'
     endif else begin
        sigma=replicate(sigma,n_dim) 
     endelse
  endif
  
  kernel = GAUSSIAN_FUNCTION(Sigma,/DOUBLE, /NORMALIZE, WIDTH=width)
  result = convol(data,kernel,_extra=extra)
  return,result
end
