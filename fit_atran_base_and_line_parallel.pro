; RK 10.10.2016
; fit a baseline with an atran model
; RK 4.1.2017
; adding idl parallel processing using parallelidl
; http://www.iluvatar.org/~dwijn/parallelidl
;RK 9.5.2017 modfied to base_and_line

function get_chunk,cube,inst,chunks
  s=size(cube)
  ch_size=float(s[2])/float(chunks)
  if inst lt chunks-1 then begin
     return,cube[*,round(inst*ch_size):round((inst+1)*ch_size)-1,*]
  endif else begin
     return,cube[*,round(inst*ch_size):*,*]     
  endelse
end

pro insert_chunk,target,source,inst,chunks
  st=size(target)
  ss=size(source)
  ch_size=float(st[2])/float(chunks)

  case st[0] of
     2: begin
        if inst lt chunks-1 then begin
           target[*,round(inst*ch_size):round((inst+1)*ch_size)-1]=source
        endif else begin
           target[*,round(inst*ch_size):*]=source
        endelse
     end
     3: begin
        if inst lt chunks-1 then begin
           target[*,round(inst*ch_size):round((inst+1)*ch_size)-1,*]=source
        endif else begin
           target[*,round(inst*ch_size):*,*]=source
        endelse
     end
  endcase
end

function fix_chunk,fix,instance,chunks
  if n_elements(fix) gt 3 $
  then return, get_chunk(fix,instance,chunks) $
  else return, fix
end

;fn_cube: filename of the cube
;out_prefix: prefix for the output files
;fit_wavel: array of wavelength interval denoting the fitting range
;Keywords:
; sub_cube: 2x2 array with indices to select subcube to be fitted - [[x0,x1],[y0,y1]]
; cube_ext: extension of fits file with fluxes to fit (defaut 3)
; err_ext: extension of fits file with flux uncertainties (defaut 4)
; scale_err: rather than using an extenstion for the uncertainties,
;            scale the flux with this value to estimate errors
; no_parallel: no parallel processing
; cpus: how many cpus to use (default 8)
; jobs_per_cpu: How often to task each cpu (default 4)

pro fit_atran_base_and_line_parallel,fn_cube,out_prefix,fit_wavel,sub_cube=sub_cube,$
                                     cube_ext=cube_ext,err_ext=err_ext,scale_err=scale_err,$
                                     no_parallel=no_parallel,cpus=cpus,jobs_per_cpu=jobs_per_cpu,$
                                     fix_alt=fix_alt,fix_pwv=fix_pwv,fix_zero=fix_zero,$
                                     fix_center=fix_center,fix_width=fix_width,$
                                     fix_resolution=fix_resolution,_extra=extra

  if n_elements(cube_ext) eq 0 then cube_ext=3
  if n_elements(err_ext) eq 0 then err_ext=4

  dummy=mrdfits(fn_cube,0,header)
  cube=mrdfits(fn_cube,cube_ext,hd_cube)
  if n_elements(sub_cube) eq 4 then begin
     cube=cube[sub_cube[0,0]:sub_cube[1,0],sub_cube[0,1]:sub_cube[1,1],*]
  endif
  check_fits,cube,header,/update

  s=size(cube)
  
  if n_elements(scale_err) eq 0 then begin
     err_cube=mrdfits(fn_cube,err_ext)
     if n_elements(sub_cube) eq 4 then begin
        err_cube=err_cube[sub_cube[0,0]:sub_cube[1,0],sub_cube[0,1]:sub_cube[1,1],*]
     endif
  endif else err_cube=scale_err*cube
  
  wave=wave_hdr(header)
  
  if ~keyword_set(cpus) then cpus=4
  if ~keyword_set(jobs_per_cpu) then jobs_per_cpu=8
  chunks = jobs_per_cpu*cpus < s[2]
  if ~keyword_set(no_parallel) then bridges = build_bridges(cpus)

  results=hash()
  fit_wavel_idx=where(fit_wavel[0] le wave and wave le fit_wavel[1])
  

  
  for instance=0,chunks-1 do begin
     print,"Chunk ",instance+1," of ", chunks
     command="fit_atran_base_and_line,cube,err_cube,header,fit_wavel_idx,"+$
             "alt_fit=alt_fit,pwv_fit=pwv_fit,cont_fit=cont_fit,zero_fit=zero_fit,"+$
             "flux_fit=flux_fit,center_fit=center_fit,width_fit=width_fit,res_fit=res_fit,"+$
             "alt_sig=alt_sig,pwv_sig=pwv_sig,cont_sig=cont_sig,zero_sig=zero_sig,"+$
             "flux_sig=flux_sig,center_sig=center_sig,width_sig=width_sig,res_sig=res_sig,"+$
             "residuals=residuals,fix_alt=fix_alt,fix_pwv=fix_pwv,fix_zero=fix_zero,"+$
             "fix_center=fix_center,fix_width=fix_width,fix_resolution=fix_resolution,_extra=extra"
     
     userdata={instance:instance,pout:ptr_new(/allocate_heap)}
     results[instance]={pout:userdata.pout}

     fix_pwv_chunk = fix_chunk(fix_pwv,instance,chunks)
     fix_zero_chunk = fix_chunk(fix_zero,instance,chunks)
     fix_center_chunk = fix_chunk(fix_center,instance,chunks)
     fix_width_chunk = fix_chunk(fix_width,instance,chunks)
     fix_resolution_chunk = fix_chunk(fix_resolution,instance,chunks)
     
     if ~keyword_set(no_parallel) then begin
        bridge = get_idle_bridge(bridges)
        bridge->setproperty, userdata=userdata
        bridge->setvar, 'cube', get_chunk(cube,instance,chunks)
        bridge->setvar, 'err_cube', get_chunk(err_cube,instance,chunks)
        bridge->setvar, 'header', header
        bridge->setvar, 'fit_wavel_idx', fit_wavel_idx
        bridge->setvar, 'fix_pwv', fix_pwv_chunk
        bridge->setvar, 'fix_zero', fix_zero_chunk
        bridge->setvar, 'fix_center', fix_center_chunk
        bridge->setvar, 'fix_width', fix_width_chunk
        bridge->setvar, 'fix_resolution', fix_resolution_chunk

        struct_pass,extra,bridge
        bridge->setproperty, callback='fit_atran_base_and_line_callback'
     
        bridge->execute, /nowait, command
     endif else begin
        fit_atran_base_and_line,get_chunk(cube,instance,chunks),$
                                get_chunk(err_cube,instance,chunks),header,fit_wavel_idx,$
                                alt_fit=alt_fit,pwv_fit=pwv_fit,cont_fit=cont_fit,zero_fit=zero_fit,$
                                flux_fit=flux_fit,center_fit=center_fit,width_fit=width_fit,res_fit=res_fit,$
                                alt_sig=alt_sig,pwv_sig=pwv_sig,cont_sig=cont_sig,zero_sig=zero_sig,$
                                flux_sig=flux_sig,center_sig=center_sig,width_sig=width_sig,res_sig=res_sig,$
                                residuals=residuals,status=status,errmsg=errmsg,$
                                fix_alt=fix_alt,fix_pwv=fix_pwv_chunk,fix_zero=fix_zero_chunk,fix_center=fix_center_chunk,$
                                fix_width=fix_width_chunk,fix_resolution=fix_resolution_chunk,_extra=extra
        if status gt 0 then begin
           *(results[instance].pout)={status:2,alt_fit:alt_fit,pwv_fit:pwv_fit,cont_fit:cont_fit,$
                                      zero_fit:zero_fit,flux_fit:flux_fit,center_fit:center_fit,$
                                      width_fit:width_fit,res_fit:res_fit,$
                                      alt_sig:alt_sig,pwv_sig:pwv_sig,cont_sig:cont_sig,$
                                      zero_sig:zero_sig,flux_sig:flux_sig,center_sig:center_sig,$
                                      width_sig:width_sig,res_sig:res_sig,residuals:residuals}
        endif else begin
           *(results[instance].pout)={status:status,instance:userdata.instance,error:errmsg}
        endelse
     endelse
  endfor
  if ~keyword_set(no_parallel) then barrier_bridges,bridges

  alt_fit=fltarr(s[1],s[2])
  pwv_fit=fltarr(s[1],s[2])
  cont_fit=fltarr(s[1],s[2])
  zero_fit=fltarr(s[1],s[2])
  flux_fit=fltarr(s[1],s[2])
  center_fit=fltarr(s[1],s[2])
  width_fit=fltarr(s[1],s[2])
  res_fit=fltarr(s[1],s[2])
  alt_sig=fltarr(s[1],s[2])
  pwv_sig=fltarr(s[1],s[2])
  cont_sig=fltarr(s[1],s[2])
  zero_sig=fltarr(s[1],s[2])
  flux_sig=fltarr(s[1],s[2])
  center_sig=fltarr(s[1],s[2])
  width_sig=fltarr(s[1],s[2])
  res_sig=fltarr(s[1],s[2])
  residuals=fltarr(s[1],s[2],s[3])

;combine results
  for i=0,chunks-1 do begin
     if (*(results[i].pout)).status eq 2 then begin
        insert_chunk,   alt_fit,(*(results[i].pout)).alt_fit,i,chunks
        insert_chunk,   pwv_fit,(*(results[i].pout)).pwv_fit,i,chunks
        insert_chunk,  cont_fit,(*(results[i].pout)).cont_fit,i,chunks
        insert_chunk,  zero_fit,(*(results[i].pout)).zero_fit,i,chunks
        insert_chunk,  flux_fit,(*(results[i].pout)).flux_fit,i,chunks
        insert_chunk,center_fit,(*(results[i].pout)).center_fit,i,chunks
        insert_chunk, width_fit,(*(results[i].pout)).width_fit,i,chunks
        insert_chunk,   res_fit,(*(results[i].pout)).res_fit,i,chunks
        insert_chunk,   alt_sig,(*(results[i].pout)).alt_sig,i,chunks
        insert_chunk,   pwv_sig,(*(results[i].pout)).pwv_sig,i,chunks
        insert_chunk,  cont_sig,(*(results[i].pout)).cont_sig,i,chunks
        insert_chunk,  zero_sig,(*(results[i].pout)).zero_sig,i,chunks
        insert_chunk,  flux_sig,(*(results[i].pout)).flux_sig,i,chunks
        insert_chunk,center_sig,(*(results[i].pout)).center_sig,i,chunks
        insert_chunk, width_sig,(*(results[i].pout)).width_sig,i,chunks
        insert_chunk,   res_sig,(*(results[i].pout)).res_sig,i,chunks
        insert_chunk,residuals,(*(results[i].pout)).residuals,i,chunks
     endif else begin
        print,"Instance ",(*(results[i].pout)).instance," returned: ",$
              (*(results[i].pout)).status," ",(*(results[i].pout)).error
     endelse
  endfor
  base_and_lines=cube-residuals

  burn_bridges,bridges

  sxaddhist,'Data fitted by fit_atran_base_and_line_parallel',header
  sxaddhist,['Uses ATRAN to model the atmosperic transmission',$
             'ATRAN was developed and kindly provided to the SOFIA program by',$
             'Steve Lord. The reference for ATRAN is:',$
             'Lord, S. D., 1992, NASA Technical Memorandum 103957.',$
             'The 2000/2001 version of the HITRAN database was used.'],$
            header,/comment 

  header_bal=header
  sxaddhist,'Contains the baseline and line',header_bal
  writefits,out_prefix+'_base_and_lines.fits',base_and_lines,header_bal

  header_residual=header
  sxaddhist,'Contains the residuals',header_residual
  writefits,out_prefix+'_residuals.fits',residuals,header_residual

  check_fits,alt_fit,header,/update

  header_alt=header
  sxaddhist,'Contains the altitude',header_alt
  sxaddpar,header_alt,'BUNIT','feet'
  writefits,out_prefix+'_alt.fits',alt_fit,header_alt
  writefits,out_prefix+'_alt.fits',alt_sig,/append

  header_pwv=header
  sxaddhist,'Contains the water varpor',header_pwv
  sxaddpar,header_pwv,'BUNIT','micron'
  writefits,out_prefix+'_pwv.fits',pwv_fit,header_pwv
  writefits,out_prefix+'_pwv.fits',pwv_sig,/append

  header_cont=header
  sxaddhist,'Contains the continuum',header_cont
  writefits,out_prefix+'_cont.fits',cont_fit,header_cont
  writefits,out_prefix+'_cont.fits',cont_sig,/append

  header_zero=header
  sxaddhist,'Contains the zero point',header_zero
  writefits,out_prefix+'_zero.fits',zero_fit,header_zero
  writefits,out_prefix+'_zero.fits',zero_sig,/append

  header_flux=header
  sxaddhist,'Contains the line flux',header_flux
  sxaddpar,header_flux,'BUNIT',strtrim(sxpar(header,'BUNIT'))+' '+strtrim(sxpar(header,'CUNIT3'))
  writefits,out_prefix+'_flux.fits',flux_fit,header_flux
  writefits,out_prefix+'_flux.fits',flux_sig,/append

  header_center=header
  sxaddhist,'Contains the line center',header_center
  sxaddpar,header_center,'BUNIT',sxpar(header,'CUNIT3')
  writefits,out_prefix+'_center.fits',center_fit,header_center
  writefits,out_prefix+'_center.fits',center_sig,/append

  header_width=header
  sxaddhist,'Contains the line width',header_width
  sxaddpar,header_width,'BUNIT',sxpar(header,'CUNIT3')
  writefits,out_prefix+'_width.fits',width_fit,header_width
  writefits,out_prefix+'_width.fits',width_sig,/append

  header_res=header
  sxaddhist,'Contains the resolution',header_res
  sxaddpar,header_res,'BUNIT','R=l/dl'
  writefits,out_prefix+'_res.fits',res_fit,header_res
  writefits,out_prefix+'_res.fits',res_sig,/append

end
