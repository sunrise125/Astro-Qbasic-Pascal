unit RIFCART;

  interface

     uses
         Graph;


     var
        XORIG, YORIG : integer;

     Function Ascissa(X:real):integer;
     Function Ordinata(Y:real):integer;
     procedure GRAFICA;

  implementation

    procedure GRAFICA;
      var
        grDriver : Integer;
        grMode   : Integer;
        ErrCode  : Integer;
      begin
        grDriver := Detect;
        InitGraph(grDriver,grMode,'c:\tp\bgi');
        ErrCode := GraphResult;
        if ErrCode <> grOk then
            begin
              WriteLn('Graphics error :',GraphErrorMsg(ErrCode));
              Halt;
            end
        else
          begin
            XORIG := getmaxx div 2;
            YORIG := getmaxy div 2;
          end
      end;


 {le funzioni Ascissa e Ordinata trasformano un punto del sistema di
  riferimento in un punto di coordinate-schermo}

    Function Ascissa(X:real):integer;
      begin
        Ascissa := round(X + XORIG)
      end;

    Function Ordinata(Y:real):integer;
      begin
        Ordinata := round(YORIG - Y)
      end;

  end.
