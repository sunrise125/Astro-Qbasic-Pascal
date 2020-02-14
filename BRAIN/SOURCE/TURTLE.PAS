unit turtle;
  interface
    uses
      Graph, Rifcart;

    
    procedure MOSTRA;
    procedure NASCONDI;
    procedure PENNASU;
    procedure PENNAGIU;
    procedure COLORE(X:word);
    procedure INITPOS(X,Y:real);
    procedure INITANG(ANG:real);
    procedure CASA;
    procedure TARTAMOD;
    procedure TARTAEND;
    procedure AVANTI(X:real);
    procedure INDIETRO(X:real);
    procedure DESTRA(ANG:real);
    procedure SINISTRA(ANG:real);

    function MAXX:integer;
    function MAXY:integer;

    function XVAL:integer;
    function YVAL:integer;
    function COLCOR:word;
    function ANGCOR:integer;

  implementation

    const
      BM:boolean = false;
      PS:boolean = false;

    var VERSO, X0, Y0 : real;
        OLDCOLOR, SIZE : word;
        P:pointer;

    procedure PENNASU;
      begin
        PS:=true;
        OLDCOLOR:=getcolor;
        setcolor(getbkcolor);
      end;

    procedure PENNAGIU;
      begin
        PS:=false;
        setcolor(OLDCOLOR);
      end;

    procedure COLORE(X:word);
      begin
        OLDCOLOR:=X;
        if PS then
          setcolor(getbkcolor)
        else
          setcolor(OLDCOLOR)
      end;

    procedure INITPOS(X,Y:real);
      begin
        if BM then
          begin
            NASCONDI;
            X0:=X;Y0:=Y;
            moveto(Ascissa(X0),Ordinata(Y0));
            MOSTRA
          end
        else
          begin
            X0:=X;Y0:=Y;
            moveto(Ascissa(X0),Ordinata(Y0));
          end
      end;

    procedure INITANG(ANG:real);
      begin
        VERSO:=ANG;
      end;

    procedure CASA;
      begin
        INITPOS(0,0);
        INITANG(0);
      end;

    procedure TARTAMOD;
      begin
        grafica;
        casa
      end;

    procedure TARTAEND;
      begin
        closegraph
      end;

    function DIREZIONE:real;
      begin
        Direzione:=VERSO*PI/180;
      end;

    procedure AVANTI(X:real);
      begin
        if BM then
          begin
            NASCONDI;
            X0:=X*cos(DIREZIONE)+X0;
            Y0:=X*sin(DIREZIONE)+Y0;
            if PS then
              moveto(Ascissa(X0),Ordinata(Y0))
            else
              lineto(Ascissa(X0),Ordinata(Y0));
            MOSTRA;
          end
        else
          begin
            X0:=X*cos(DIREZIONE)+X0;
            Y0:=X*sin(DIREZIONE)+Y0;
            if PS then
              moveto(Ascissa(X0),Ordinata(Y0))
            else
              lineto(Ascissa(X0),Ordinata(Y0));
          end
      end;

    procedure INDIETRO(X:real);
      begin
        if BM then
          begin
            NASCONDI;
            X0:=X*cos(DIREZIONE+PI)+X0;
            Y0:=X*sin(DIREZIONE+PI)+Y0;
            if PS then
              moveto(Ascissa(X0),Ordinata(Y0))
            else
              lineto(Ascissa(X0),Ordinata(Y0));
            MOSTRA;
          end
        else
          begin
            X0:=X*cos(DIREZIONE+PI)+X0;
            Y0:=X*sin(DIREZIONE+PI)+Y0;
            if PS then
              moveto(Ascissa(X0),Ordinata(Y0))
            else
              lineto(Ascissa(X0),Ordinata(Y0));
          end
      end;

    procedure DESTRA(ANG:real);
      begin
        if BM then
          begin
            NASCONDI;
            VERSO:=VERSO-ANG;
            MOSTRA;
          end
        else
          VERSO:=VERSO-ANG;
      end;

    procedure SINISTRA(ANG:real);
      begin
        if BM then
          begin
            NASCONDI;
            VERSO:=VERSO+ANG;
            MOSTRA;
          end
        else
          VERSO:=VERSO+ANG;
      end;

    procedure MOSTRA;

      var I:byte;
          VOLD, XOLD, YOLD:real;
          OLDPS:boolean;
          X1, X2, Y1, Y2, XMAX, YMAX:integer;

      function INTERS(XA,YA,XB,YB:integer):boolean;
      {se il rettangolo di vertici opposti A(XA,YA),B(XB,YB) non
       interseca il rettangolo dello schermo grafico attivo la
       funzione restituisce false, vero altrimenti}
        begin
          if ((not(XA in [0..XMAX])) and (not(XB in [0..XMAX]))) or
             ((not(YA in [0..YMAX])) and (not(YB in [0..YMAX])))
             then INTERS:=false
             else INTERS:=true
        end;

      begin
        {memorizzo un quadrato contenente il vecchio disegno}
        XOLD:=X0;YOLD:=Y0;
        VOLD:=VERSO;
        XMAX:=Getmaxx;YMAX:=Getmaxy;
        X1:=Ascissa(XOLD-25);Y1:=Ordinata(YOLD+25);
        X2:=Ascissa(XOLD+25);Y2:=Ordinata(YOLD-25);
        {se il rettangolo Š completamente fuori dallo schermo
         esco dalla procedura non essendo comunque possibile
         visualizzare la tartaruga}
        if INTERS(X1,Y1,X2,Y2) then
          begin
            if X1<0
              then X1:=0;
            if X2>XMAX
              then X2:=XMAX;
            if Y1<0
              then Y1:=0;
            if Y2>YMAX
              then Y2:=YMAX;
            SIZE:=ImageSize(X1,Y1,X2,Y2);
            GetMem(P, SIZE); { Preleva  memoria dall' heap }
            GetImage(X1,Y1,X2,Y2,P^);
            {disegno la tartaruga}
            OLDPS:=PS;
            if not(PS) then pennasu;
            BM:=false;
            avanti(8);
            sinistra(150);
            pennagiu;
            AVANTI(15);
            SINISTRA(120);
            pennasu;
            AVANTI(15);
            pennagiu;
            SINISTRA(120);
            AVANTI(15);
            if OLDPS then pennasu;
            X0:=XOLD;Y0:=YOLD;
            moveto(Ascissa(X0),Ordinata(Y0));
            Initang(VOLD);
          end
        else SIZE:=0;
        BM:=true;
      end;

    procedure NASCONDI;
      var I:byte;
          VOLD, XOLD, YOLD:real;
          COLD :word;
          X1,Y1:integer;
      begin
        BM:=false;
        if SIZE>0 then
          begin
            X1:=Ascissa(X0-25);Y1:=Ordinata(Y0+25);
            PutImage(X1,Y1,P^,Normalput);
            freemem(P,SIZE);
          end
      end;

    {funzioni di osservazione}

    function MAXX:integer;
      begin
        MAXX:=XORIG
      end;

    function MAXY:integer;
      begin
        MAXY:=YORIG
      end;

    function XVAL:integer;
      begin
        XVAL:=round(X0)
      end;

    function YVAL:integer;
      begin
        YVAL:=round(Y0)
      end;

    function COLCOR:word;
      begin
        COLCOR:=OLDCOLOR
      end;

    function ANGCOR:integer;
      begin
        ANGCOR:=round(VERSO)
      end;

  end.
