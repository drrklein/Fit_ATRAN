; RK 10.10.2016
; fit a baseline with an atran model
; RK 4.1.2017
; adding idl parallel processing using parallelidl
; http://www.iluvatar.org/~dwijn/parallelidl
; RK 15.5.2017
; adding a continuum slope
; RK 12.9.2017
; removing emission which is no longer supported by fit_atran_baseline

function get_chunk,cube,inst,chunks
  s=size(cube)
  ch_size=s[2]/chunks
  if inst lt chunks-1 then begin
     return,cube[*,inst*ch_size:(inst+1)*ch_size-1,*]
  endif else begin
     return,cube[*,inst*ch_size:*,*]     
  endelse
end

pro insert_chunk,target,source,inst,chunks
  st=size(target)
  ss=size(source)
  ch_size=st[2]/chunks

  case st[0] of
     2: begin
        if inst lt chunks-1 then begin
           target[*,inst*ch_size:(inst+1)*ch_size-1]=source
        endif else begin
           target[*,inst*ch_size:*]=source
        endelse
     end
     3: begin
        if inst lt chunks-1 then begin
           target[*,inst*ch_size:(inst+1)*ch_size-1,*]=source
        endif else begin
           target[*,inst*ch_size:*,*]=source
        endelse
     end
  endcase
end

;fn_cube: filename of the cube
;out_prefix: prefix for the output files
;base_wavel: [2,N] array of wavelength intervals denoting the
;      baseline fitting range
pro fit_atran_baseline_parallel,fn_cube,out_prefix,base_wavel,cube_ext=cube_ext,err_ext=err_ext,$
              pwv_fit=pwv_fit,cont_fit=cont_fit,zero_fit=zero_fit,slope_fit=slope_fit,$
              bl_subtracted=bl_subtracted,baselines=baselines,za=za,alt=alt,$
              clip=clip,no_zero=no_zero,fix_pwv=fix_pwv,fix_slope=fix_slope,$
              no_parallel=no_parallel,cpus=cpus,jobs_per_cpu=jobs_per_cpu
                                
                                
  if n_elements(cube_ext) eq 0 then cube_ext=3
  if n_elements(err_ext) eq 0 then err_ext=4
  if ~keyword_set(clip) then clip=0.
  clip_str=string(clip)
  if keyword_set(no_zero) then no_zero=',/no_zero' else no_zero=''
  if ~keyword_set(fix_pwv) then fix_pwv=0.
  fix_pwv_str=string(fix_pwv)
  case n_elements(fix_slope) of
     0:fix_slope_str=''
     1:fix_slope_str=string(fix_slope,format='(",fix_slope=",E12.5)')
     2:fix_slope_str=string(fix_slope,format='(",fix_slope=[",E12.5,",",E12.5,"]")')
  endcase
  
  dummy=mrdfits(fn_cube,0,header)
  cube=mrdfits(fn_cube,cube_ext,hd_cube)
  check_fits,cube,header,/update
  err_cube=mrdfits(fn_cube,err_ext)

  wave=wave_hdr(header)
  base_idx=where(base_wavel[0,0] le wave and wave le base_wavel[1,0])
  for i=1,n_elements(base_wavel)/2-1 do begin
     base_idx=[base_idx,where(base_wavel[0,i] le wave and wave le base_wavel[1,i])]
  endfor
  
  get_za_alt,header,za,alt; just to pass it back

  if ~keyword_set(cpus) then cpus=4
  if ~keyword_set(jobs_per_cpu) then jobs_per_cpu=8
  chunks=jobs_per_cpu*cpus
  if ~keyword_set(no_parallel) then bridges = build_bridges(cpus)
  results=hash()
  
  s=size(cube)
  
  for instance=0,chunks-1 do begin
     print,"Chunk ",instance+1," of ", chunks
     command="fit_atran_baseline,cube,err_cube,header,base_idx,"+ $
             "pwv_fit=pwv_fit,cont_fit=cont_fit,zero_fit=zero_fit,slope_fit=slope_fit,"+ $
             "bl_subtracted=bl_subtracted,/quiet,clip="+clip_str+ $
             ",fix_pwv="+fix_pwv_str+no_zero+fix_slope_str
     
     if ~keyword_set(no_parallel) then begin
        bridge = get_idle_bridge(bridges)
        userdata={instance:instance,pout:ptr_new(/allocate_heap)}
        results[instance]={pout:userdata.pout}
        bridge->setproperty, userdata=userdata
     
        bridge->setvar, 'cube', get_chunk(cube,instance,chunks)
        bridge->setvar, 'err_cube', get_chunk(err_cube,instance,chunks)
        bridge->setvar, 'header', header
        bridge->setvar, 'base_idx', base_idx
        bridge->setproperty, callback='fit_atran_baseline_callback'
     
        bridge->execute, /nowait, command
     endif else begin
        result=execute(command)
     endelse
  endfor
  if ~keyword_set(no_parallel) then barrier_bridges,bridges
  
  pwv_fit=fltarr(s[1],s[2])
  cont_fit=fltarr(s[1],s[2])
  zero_fit=fltarr(s[1],s[2])
  slope_fit=fltarr(s[1],s[2])
  bl_subtracted=fltarr(s[1],s[2],s[3])
  
;combine results
  for i=0,chunks-1 do begin
     if (*(results[i].pout)).status eq 2 then begin
        insert_chunk,  pwv_fit,(*(results[i].pout)).pwv_fit,i,chunks
        insert_chunk, cont_fit,(*(results[i].pout)).cont_fit,i,chunks
        insert_chunk, zero_fit,(*(results[i].pout)).zero_fit,i,chunks
        insert_chunk,slope_fit,(*(results[i].pout)).slope_fit,i,chunks
        insert_chunk,bl_subtracted,(*(results[i].pout)).bl_subtracted,i,chunks
     endif else begin
        print,"Instance ",(*(results[i].pout)).instance," returned: ",$
              (*(results[i].pout)).status," ",(*(results[i].pout)).error
     endelse
  endfor
  baselines=cube-bl_subtracted
  
  if ~keyword_set(no_parallel) then burn_bridges,bridges

  sxaddhist,'Created by fit_atran_baseline_parallel',header

  header_bl=header
  sxaddhist,'Contains the baselines',header_bl
  writefits,out_prefix+'_baselines.fits',baselines,header_bl
  
  header_blsbt=header
  sxaddhist,'Contains the baseline subtracted cube',header_blsbt
  writefits,out_prefix+'_bl_subtracted.fits',bl_subtracted,header_blsbt
  
  check_fits,pwv_fit,header,/update

  header_pwv=header
  sxaddhist,'Contains the precipitbale water vapor',header_pwv
  writefits,out_prefix+'_pwv.fits',pwv_fit,header_pwv
  
  header_c=header
  sxaddhist,'Contains the continuum',header_c
  writefits,out_prefix+'_cont.fits',cont_fit,header_c
  
  header_z=header
  sxaddhist,'Contains the zero point',header_z
  writefits,out_prefix+'_zero.fits',zero_fit,header_z
  
  header_s=header
  sxaddhist,'Contains the slope',header_s
  writefits,out_prefix+'_slope.fits',slope_fit,header_s

end
