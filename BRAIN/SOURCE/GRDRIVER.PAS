Unit GrDriver;
INTERFACE
USES GRAPH;
IMPLEMENTATION

PROCEDURE AttDriverProc; External ;{$L Att.Obj}
PROCEDURE CGADriverProc; External ;{$L CGA.obj}
PROCEDURE EGAVGADriverProc; External ;{$L EGAVGA.OBJ}
PROCEDURE HERCDriverProc; External ; {$L HERC.OBJ}

PROCEDURE ReportError(S:String);
BEGIN
 Writeln;
 writeln(S,':',GraphErrorMsg(GraphResult));
 Halt(1);
END;

BEGIN
 if RegisterBGIDriver(@ATTDriverProc)<0
  then ReportError('AT&T');
 if RegisterBGIDriver(@CGADriverProc)<0
  then ReportError('CGA');
 if RegisterBGIDriver(@EGAVGADriverProc)<0
  then ReportError('AT&T');
 if RegisterBGIDriver(@ATTDriverProc)<0
  then ReportError('Ega-Vga');
 if RegisterBGIDriver(@HercDriverProc)<0
  then ReportError('HERC');
END.
