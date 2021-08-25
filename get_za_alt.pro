;RK 13.4.2017
;get zenith angle and altitude from SOFIA-fits headers e.g. for atran models
pro get_za_alt,header,za,alt,start=start,stop=stop
  
zasta = float(fxpar(header, 'ZA_START'))
zaend = float(fxpar(header, 'ZA_END'))
altsta = float(fxpar(header, 'ALTI_STA'))
altend = float(fxpar(header, 'ALTI_END'))

case 1 of
   keyword_set(start) and ~keyword_set(stop): begin
      za = zasta
      alt = altsta
   end
   ~keyword_set(start) and keyword_set(stop): begin
      za = zaend 
      alt = altend 
   end
   keyword_set(start) and keyword_set(stop): begin
      za=[zasta,zaend]
      alt=[altsta,altend]
   end
   ~keyword_set(start) and ~keyword_set(stop):begin
      za = (zasta + zaend) / 2.0
      alt = (altsta + altend) / 2.0
   end
   else: message,"You should't come here"
endcase

end
