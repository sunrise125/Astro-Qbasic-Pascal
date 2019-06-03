' Calcolo del Tempo siderale a Greenwich
' =======================================  (Ott.04)
  COLOR 14, 1: CLS
  DEFDBL A-Z
  pi = 4 * ATN(1): rad = pi / 180
'------------------------------------------------------------------------------------------------------------------------------------
COLOR 14, 4
  LOCATE 1, 4: PRINT "здддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддд©"
  LOCATE 2, 4: PRINT "Ё                                                                        Ё"
  LOCATE 3, 4: PRINT "юдддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддддды"
COLOR 14, 5
  LOCATE 2, 6: PRINT "               Calcolo del TEMPO SIDERALE a GREENWICH                 "
COLOR 14, 1
  LOCATE 5, 17: PRINT "здддддддддддд Data di riferimento дддддддддддд©"
  LOCATE 6, 17: PRINT "Ё                                             Ё"
  LOCATE 7, 17: PRINT "Ё                                             Ё"
  LOCATE 8, 17: PRINT "юддддддддддддддддддддддддддддддддддддддддддддды"

 DO: LOCATE 6, 20: INPUT "Giorno: ", gio
      LOOP UNTIL gio <= 31 AND gio > 0
 DO: LOCATE 6, 37: INPUT "Mese: ", mes
      LOOP UNTIL mes <= 12 AND mes > 0
     LOCATE 6, 51: INPUT "Anno: ", ann

 DO: LOCATE 7, 21: INPUT "Ora (Greenw.o Solare o Legale) (g/s/l) ", wo$
      LOOP WHILE wo$ <> "s" AND wo$ <> "l" AND wo$ <> "g"
        IF wo$ = "g" THEN cc = 0: 'Fuso orario 0
        IF wo$ = "s" THEN cc = 1: 'Fuso orario 1
        IF wo$ = "l" THEN cc = 2: 'Fuso orario 2
'---------------------------------------------------------------------------

  LOCATE 9, 32: PRINT "Fuso Orario: "; cc

  LOCATE 11, 24: PRINT "здддддддд    Orario   ддддддддд©"
  LOCATE 12, 24: PRINT "Ё Ore:     Min:      Sec:      Ё"
  LOCATE 13, 24: PRINT "юдддддддддддддддддддддддддддддды"


 'Calcolo di JD
      u7 = 2 - INT(ann / 100) + INT(INT(ann / 100) / 4)
      IF mes < 3 THEN
          ann = ann - 1: mes = mes + 12
          END IF
      jd = INT(365.25 * ann) + INT(30.6001 * (mes + 1)) + gio + 1720994.5#
      IF jd >= 2299160.5# THEN jd = jd + u7
      jq = jd - 2451545: t = jq / 36525

' -------------------------
 DO:  LOCATE 12, 31: INPUT "", ore
       LOOP UNTIL ore <= 23 AND ore >= 0
 DO:  LOCATE 12, 40: INPUT "", min
       LOOP UNTIL min <= 59 AND min >= 0
 DO:  LOCATE 12, 50: INPUT "", sec
       LOOP UNTIL sec <= 59 AND sec >= 0

  GOSUB subr

  LOCATE 16, 12:
  PRINT " Tempo Siderale a Greenwich:";
  PRINT USING " ##.###### ore"; tsidd;
  delgr# = tsidd
  GOSUB GradAsess
  PRINT USING "= ##h ##m ##.##s"; deltagra#; deltamin#; deltasec#
  PRINT
  PRINT USING "                                     = ###.##### gradi"; tsidd * 15
 END

'------------- SUBROUTINES ------------------------------------------------
subr:
 'Calcolo del T-sid a Greenwich all'ora fissata -
      wn1 = (ore - cc) + min / 60 + sec / 3600: ' tempo universale TU
     mat1 = .279057273# + 100.002139027# * t + 1.0776E-06 * t ^ 2
     mat2 = mat1 - INT(mat1)
   tsidd = 1.00273790935# * wn1 + mat2 * 24
   tsidd = tsidd - INT(tsidd / 24) * 24: 'in ORE

    ts = tsidd * 15: 'in GRADI

RETURN

GradAsess:
' Conversione di GRADI in GRADI,PRIMI,SEC e GRADI,PRIMI,Fraz.PRIMI
'
 deltagra# = INT(ABS(delgr#))
 deltafraz# = (ABS(delgr#) - deltagra#) * 60#
 deltamin# = INT(deltafraz#)
 deltasec# = (deltafraz# - deltamin#) * 60#
 deltagra# = deltagra# * SGN(delgr#)
'
 RETURN

' ******** Fine Programma: TSID.BAS ************ (Ott.04) ****************

