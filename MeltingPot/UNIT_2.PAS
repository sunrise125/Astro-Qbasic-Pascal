 Unit UNIT_2;    { 2^ Unit… per Codice SWIFT
***********************************************************************
c                         ORBEL_XV2EL.F
***********************************************************************
*     PURPOSE:  Given the cartesian position and velocity of an orbit,
*       compute the osculating orbital elements.
*
C       input:
c            x,y,z    ==>  position of object (real scalars)
c            vx,vy,vz ==>  velocity of object (real scalars)
c            gmsum       ==> G*(M1+M2) (real scalar)
c
c       Output:
c      ialpha   ==> conic section type ( see PURPOSE, integer scalar)
C      a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
c            omega    ==> argument of perihelion (real scalar)
C            capm     ==> mean anomoly(real scalar)
c
*     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech." 
*     REMARKS:  If the inclination INC is less than TINY, we
*       arbitrarily choose the longitude of the ascending node LGNODE
*       to be 0.0 (so the ascending node is then along the X axis).  If 
*       the  eccentricity E is less than SQRT(TINY), we arbitrarily
*       choose the argument of perihelion to be 0.
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  May 8,1992.
*     REVISIONS:
***********************************************************************
  subroutine orbel_xv2el(x,y,z,vx,vy,vz,gmsum,
     &     ialpha,a,e,inc,capom,omega,capm)
}

 Interface

 Uses Crt;
 Const  TINY=4e-15; rad=pi/180;
 Var
      {Solo Internals e Outputs; gli Inputs definiti nel Progr. principale
       -------------------------------------------------------------------}
     ialpha:   Integer;
     a,e,inc,capom,omega,capm,
     hx,hy,hz,h2,h,r,v2,v,vdotr,energy,fac,face,cape,capf,tmpf,
     cw,sw,w,u
                                 : Real;

 Procedure ORBEL_XV2EL(x,y,z,vx,vy,vz,gmsum: Real);
      Function ArcoCOS(valore:Real): Real;
      Function ArcoTANG(num,den:Real): Real;

 Implementation

      Function ArcoCOS(valore:Real): Real;
   Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(abs(1-valore*valore)));
  End;

      Function ArcoTANG(num,den:Real): Real;
    Var atan: Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

      Function COSH(valore:Real): Real;
   Begin
     COSH:=(exp(valore)+exp(-valore))/2;
  End;

      Function SINH(valore:Real): Real;
   Begin
     SINH:=(exp(valore)-exp(-valore))/2;
  End;

 Procedure ORBEL_XV2EL(x,y,z,vx,vy,vz,gmsum: Real);
   Begin

 { Compute the angular momentum H, and thereby the inclination INC. }

       hx:= y*vz - z*vy;
       hy:= z*vx - x*vz;
       hz:= x*vy - y*vx;
       h2:= hx*hx + hy*hy +hz*hz;
        h:= sqrt(h2);
      inc:= ArcoCOS(hz/h);

 { Compute longitude of ascending node CAPOM and the argument of latitude u. }
       fac:= sqrt(hx*hx + hy*hy)/h;

       if fac< TINY then begin
         capom:= 0;
         u:= ArcoTANG(y,x); end;

       if abs(inc - PI) < 10*TINY then u:= -u
       else begin
         capom:= arcoTANG(hx,-hy);
         u:= arcoTANG ( z/sin(inc) , x*cos(capom) + y*sin(capom));
       end;

       if capom < 0  then capom:= capom + 2*PI;
       if u < 0  then u:= u + 2*PI;


 { Compute the radius R and velocity squared V2, and the dot
     product RDOTV, the energy per unit mass ENERGY . }

       r:= sqrt(x*x + y*y + z*z);
       v2:= vx*vx + vy*vy + vz*vz;
       v:= sqrt(v2);
       vdotr:= x*vx + y*vy + z*vz;
       energy:= 0.5*v2 - gmsum/r;

{  Determine type of conic section and label it via IALPHA }

       if abs(energy*r/gmsum) < sqrt(TINY) then ialpha:= 0
       else begin
          if energy < 0 then ialpha:= -1;
          if energy > 0 then ialpha:= +1;
       end;

 { Depending on the conic type, determine the remaining elements

        ELLIPSE :  }

 if ialpha = -1 then begin
            a:= -0.5*gmsum/energy;
            fac:= 1 - h2/(gmsum*a);
   if  fac > TINY  then begin
                 e:= sqrt ( fac );
                 face:=(a-r)/(a*e);
             if face > 1 then cape:= 0
                 else if face > -1 then cape:= arcoCOS( face )
                     else begin
                       cape:= PI;
             end;
                     if vdotr < 0  then cape:= 2*pi - cape;
                     cw:= (cos( cape) -e)/(1 - e*cos(cape));
                     sw:= sqrt(1-e*e)*sin(cape)/(1 - e*cos(cape));
                     w:= arcoTANG(sw,cw);
                     if w < 0 then w:= w + 2*PI
        end
     else begin
          e:= 0;
          w:= u;
          cape:= u;
   end;
           capm:= cape - e*sin (cape);
           omega:= u - w;
           if omega < 0 then omega:= omega + 2*PI;
           omega:= omega - int(omega/(2*PI))*2*PI;
 end;

{   HYPERBOLA  }
 if ialpha = +1 then begin
           a:= +0.5*gmsum/energy;
            fac:= h2/(gmsum*a);

   if fac > TINY then begin
            e:= sqrt ( 1 + fac );
           tmpf:= (a+r)/(a*e);
           capf:= ln(tmpf + sqrt(tmpf*tmpf -1));
           if  vdotr < 0 then capf:= - capf;
           cw:= (e - COSH(capf))/(e*COSH(capf) - 1 );
           sw:= sqrt(e*e - 1)*SINH(capf)/(e*COSH(capf) - 1 );
           w:= arcoTANG(sw,cw);
           if w < 0 then w:= w + 2*PI; end
         else begin
     { we only get here if a hyperbola is essentially a parabola
       so we calculate e and w accordingly to avoid singularities }
           e:= 1;
           tmpf:= 0.5*h2/gmsum;
           w:= arcoCOS(2*tmpf/r -1);
           if vdotr < 0 then w:= 2*PI - w;
           tmpf:= (a+r)/(a*e);
           capf:= ln(tmpf + sqrt(tmpf*tmpf -1));
   end;
         capm:= e * SINH(capf) - capf;
         omega:= u - w;
         if omega < 0 then omega:= omega + 2*PI;
         omega:= omega - int(omega/(2*PI))*2*PI;
 end;

{  PARABOLA : ( NOTE - in this case we use "a" to mean pericentric distance) }
       if ialpha = 0 then begin
         a:=  0.5*h2/gmsum;
         e:= 1;
         w:= arcoCOS(2*a/r -1);
         if  vdotr < 0 then w:= 2*PI - w;
         tmpf:= sin(0.5 * w)/cos(0.5 * w);
         capm:= tmpf* (1 + tmpf*tmpf/3);
         omega:= u - w;
         if omega < 0 then omega:= omega + 2*PI;
         omega:= omega - int(omega/(2*PI))*2*PI;
       end;

 end;  {Chiude ORBEL_XV2EL}

 END.
