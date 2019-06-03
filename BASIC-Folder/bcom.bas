' https://en.wikibooks.org/wiki/QBasic/Appendix
' ---------------------------------------------
COLOR 14, 1
OPEN "lista.txt" FOR INPUT AS #1
kount = 0
  DO WHILE NOT EOF(1)
    kount = kount + 1
      LINE INPUT #1, REC$  'Read entries from the file.
      PRINT kount,
      SHELL "mpv.exe REC$"           'Print the entries on the screen.
      LOOP
CLOSE #1


