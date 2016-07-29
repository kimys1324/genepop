program  linkdos;

   { calculates linkage desequilibria between couple of loci }
   {   in several subdivided populations, and the OHTA's     }
   {   components of variance of the desequilibrium within   }
   {      and between subpopulation                          }

{$E+}
{$N+}

uses  Dos,Crt;

const     npopm=40; nlocm=20; alm=99;        {wich may be modified}

type      ptald    =^matald;
          ptdl     =^matdl;
          ptdes    =^matdes;
          ptstat   =^matstat;
          ptst     =^matst;

          matpop   = array [1..npopm] of integer;
          matal    = array [1..npopm,1..nlocm,1..alm] of single;
          matald   = array [1..alm,1..nlocm] of integer;
          matdl    = array [1..alm,1..alm] of single;
          matind   = array [1..npopm,1..nlocm] of integer;
          matnumal = array [1..nlocm] of integer;

          matdes   = array [1..alm,1..alm] of single;
          matstat  = array [1..npopm,1..nlocm] of single;
          matst    = array [1..nlocm,1..alm] of single;


var       popn,nmax                              : matpop;
          alef                                   : matal;
          alefld,homold,sumalf,homosm            : ptald;
          sumld,totsum                           : ptdl;
          ind                                    : matind;
          npopct,indld,fixchk,nodat,sumind,negfre : integer;
          nall,totn,ndr1,ndr2                    : matnumal;

          burrow,corrad                          : ptdes;
          bursq,bursd,chisq,z,prob               : single;
          fone,fthree,fseven,seven               : single;
          chisc,rcprob                           : matdl;
          wta                                    : ptst;
          apop                                   : array [1..npopm] of string[10];
          aloc                                   :array [1..nlocm] of string[5];
          input,outdes,outohta,outfreq           : string[15];
          npop,nloc                              : integer;
          signi                                  : real;

          file1,file2,file3,filet                : text;
          ldopt,ohtaopt,freqopt,ldtotopt         : byte;
          response                               : char;
          pop,loc,loco,cc                        : integer;
          allchk,df                              : integer;



function power(x,n  :real):real;

var  signe    :integer;
     parit    :integer;

 begin
   if x=0 then power:=1
   else power:=exp(n*Ln(abs(x)));
 end;



(**********************************************************************)
(**************Integrales du chi2 et de la loi normale*****************)
(**********************************************************************)


function gauss (x  : real):real;

const    zero=0.000001;
var      z,y,w        :double;

begin
  if (not (x=zero)) then begin
     y:=abs(x)/2;
     if (y>=3) then z:=1
     else if (y>1) then begin
          y:=y-2;
          z:=-0.000045155659*y+0.000152529290;
          z:=z*y-0.000019538132;
          z:=z*y-0.000676904986;
          z:=z*y+0.001390604284;
          z:=z*y-0.000794620820;
          z:=z*y-0.002034254874;
          z:=z*y+0.006549791214;
          z:=z*y-0.010557625006;
          z:=z*y+0.011630447319;
          z:=z*y-0.009279453341;
          z:=(z*y+0.005353579108)*y-0.002141268741;
          z:=(z*y+0.000535310849)*y+0.999936657524;

{          z:=(((((((((((((-0.045155659*y
              +0.152529290)*y-0.019538132)*y
              -0.676904986)*y+1.390604284)*y
              -0.794620820)*y-2.034254874)*y
              +6.549791214)*y-10.557625006)*y
              +11.630447319)*y-9.279453341)*y
              +5.353579108)*y-2.141268741)*y
              +0.535310849)*y+999.936657524;
          z:=z/1000;}
    end
    else begin
         w:=sqr(y);
         z:=0.000124818987*w-0.001075204047;
         z:=z*w+0.005198775019;
         z:=(z*w-0.019198292004)*w+0.059054035642;
         z:=(z*w-0.151968751364)*w+0.319152932694;
         z:=((z*w-0.531923007300)*w+0.797884560593)*y*2;
{         z:=((((((((0.000124818987*w
             -0.001075204047)*w+0.005198775019)*w
             -0.019198292004)*w+0.059054035642)*w
             -0.151968751364)*w+0.319152932694)*w
             -0.531923007300)*w+0.797884560593)*y*2;}
    end;
  end
  else z:=0;

  if (x>zero) then gauss:=(1+z)/2
              else gauss:=(1-z)/2;

end;    {gauss}


function  chiprb(chisq : real;df :integer)  :single;

const   zero=0.000001;
var even,bigx      : boolean;
    check          : integer;
    a,s,z,e,c      : real;
    chisqb,y       : real;

begin
  if (chisq<=zero) then chiprb:=0
  else begin
    if df>30 then chiprb:=1-gauss(-sqrt(4.5*df)*(power(chisq/df,0.333333333)
                                     +2/(9*df)-1))
    else begin
      a:=0.5*chisq;
      chiprb:=0;
      check:=2*(df div 2);
        if (check=df) then even:=true
                      else even:=false;
        if chisq>=200 then bigx:=true
                      else bigx:=false;
      if df<=2 then begin
        if (bigx=true) then chiprb:=1
        else begin
          if (even=true) then chiprb:=1-exp(-a)
                         else chiprb:=1-2*gauss(-sqrt(chisq));
        end;
      end
      else begin {df>2}
         chisqb:=0.5*(df-1);
         if (bigx=true) then  begin
            if (even=true) then begin
              s:=0; z:=0 ; e:=0;
            end
            else begin
              s:=2*gauss(-sqrt(chisq));
              z:=-0.5;
              e:=0.57236494295   {ln(sqrt(pi))};
            end;
            if z<=0 then z:=0.000001;
            c:=ln(z)+e;
            z:=z+1;
            while (z<chisqb) do begin
               e:=ln(z+e);
               s:=exp(c*z-a-e)+s;
               z:=z+1
            end;
            chiprb:=1-s
         end  {bigx=true}
         else begin {bigx=false}
           y:=exp(-a);
           if (even=true) then begin
             z:=0;
             s:=y;
             e:=1
           end
           else begin
             s:=2*gauss(-sqrt(chisq));
             z:=-0.5;
             e:=0.564189583548/sqrt(a) {1/sqrt(pi)};
           end;
           c:=0;
           z:=z+1;
           while (z<chisqb) do begin
            e:=e*a/z;
            c:=c+e;
            z:=z+1
           end;
           chiprb:=1-(c*y+s);
         end; {bigx=false}
       end;{df>2}
   end;{df<30}
  end; {chisq>zero}
end;


(**********************************************************************)
(**************Initialisation des matrices*****************************)
(**********************************************************************)

Procedure init;

var all                    :integer;
    popi,loci              :integer;


begin
   for popi:=1 to npopm do begin
      popn[popi]:=0;
      nmax[popi]:=0;
     for loci:=1 to nlocm do begin
        ind[popi,loci]:=0;
       for all:=1 to alm do begin
          alef[popi,loci,all]:=0;
          wta^[loci,all]:=0;
       end;
     end;
   end;
   for loci:=1 to nlocm do begin
        nall[loci]:=0;
        totn[loci]:=0;
   end;

end;  {init}





(**********************************************************************)
(***********  Initialisation des variables et des matrices  ***********)
(***********    pour chaque couple de locus  **************************)
(**********************************************************************)


procedure initcouple (loc1,loc2  :integer);

var  al1,al2    :integer;

begin
         npopct:=npop;
         indld:=0;
         fixchk:=0;
         nodat:=0;
         sumind:=0;
         negfre:=0;
         bursq :=0;
         bursd :=0;
         chisq :=0;
         z :=0;
         prob:=0;
         fone:=0;
         fthree:=0;
         fseven:=0;
         seven :=0;
        for al1:=1 to alm do
          for al2:=1 to alm do begin
            sumld^[al1,al2]:=0;
            totsum^[al1,al2]:=0;
            burrow^[al1,al2]:=0;
            corrad^[al1,al2]:=0;
            alefld^[al1,loc2]:=0;
            alefld^[al2,loc1]:=0;
            homold^[al1,loc2]:=0;
            homold^[al2,loc1]:=0;
            sumalf^[al1,loc2]:=0;
            sumalf^[al2,loc1]:=0;
            homosm^[al1,loc2]:=0;
            homosm^[al2,loc1]:=0;
            chisc[al1,al2]:=0;
            rcprob[al1,al2]:=0;
          end;

end;   {initcouples}



(**********************************************************************)
(****Lecture du fichier de donn‚es                           **********)
(**********************************************************************)

procedure lecture;

var  indi,i,j,loc           :integer;


begin


         for indi:=1 to nmax[pop] do
         begin
             for loc:=1 to nloc do
             begin
              read (file1,ndr1[loc],ndr2[loc]);
               for i:=1 to nall[loc] do
               begin
                 if ndr1[loc]<>0 then begin
                   if ndr1[loc]=i then begin
                     for j:=1 to nall[loc] do
                     begin
                       if ndr2[loc]<>0 then begin
                          if ndr2[loc]=j then begin
                           if i<>j then begin
                               alef[pop,loc,i]:=1+ alef[pop,loc,i];
                               alef[pop,loc,j]:=1+ alef[pop,loc,j];
                               ind[pop,loc]:=1+ind[pop,loc];
                              end
                             else if (i=j) then begin
                               alef[pop,loc,i]:=2+ alef[pop,loc,i];
                               ind[pop,loc]:=1+ind[pop,loc];
                           end;  {cas o— i<>j ou i=j}
                          end; {condition ndr2[loc]=j}
                       end; {ndr2[loc]<>0}

                     end; {boucle sur ndr2[loc]}
                   end; {condition ndr1[loc]=i}
                 end; {condition ndr1[loc]<>0}
               end; {boucle sur ndr1[loc]}
             end; {boucle sur les locus}
         end;     {boucle sur les individus intra sous-pop}

 end;   {lecture}




(**********************************************************************)
(****Remplissage des matrices                                **********)
(**********************************************************************)

procedure remp (lc1,lc2  :integer);


var  i,j,k,l             :integer;
     indi                :integer;


   {Comptage des diff‚rents types de gamŠtes par couple de locus, etc...}

begin

 for indi:=1 to nmax[pop] do begin

  for i:=1 to nloc do read (file1,ndr1[i],ndr2[i]);
  readln (file1);

  for i:=1 to nall[lc1] do   {compteur 1er allŠle 1er locus}
   begin
    if ndr1[lc1]<>0 then begin
      if ndr1[lc1]=i then begin
        for j:=1 to nall[lc1] do  {compteur 2Šme allŠle  1er locus}
        begin
          if ndr2[lc1]<>0 then begin
            if ndr2[lc1]=j then begin

               for k:=1 to nall[lc2] do  {compteur 1er allŠle 2Šme locus}
                begin
                  if ndr1[lc2]<>0 then begin
                    if ndr1[lc2]=k then begin
                      for l:=1 to nall[lc2] do {compteur 2Šme allŠle 2Šme loc}
                      begin
                        if ndr2[lc2]<>0 then begin
                          if ndr2[lc2]=l then begin

{cas o— i=j et k=l}         if ((i=j) and (k=l)) then begin

                              sumld^[i,k]:=2 + sumld^[i,k];
                              alefld^[i,lc2]:=2 + alefld^[i,lc2];
                              alefld^[k,lc1]:=2 + alefld^[k,lc1];
                              homold^[i,lc2]:=1 + homold^[i,lc2];
                              homold^[k,lc1]:=1 + homold^[k,lc1];
                              indld:= 1 + indld;

                              totsum^[i,k]:=2 +totsum^[i,k];
                              sumalf^[i,lc2]:=2 + sumalf^[i,lc2];
                              sumalf^[k,lc1]:=2 + sumalf^[k,lc1];
                              homosm^[i,lc2]:=1 + homosm^[i,lc2];
                              homosm^[k,lc1]:=1 + homosm^[k,lc1];
                              sumind := 1 + sumind;

                            end;

{cas o— i=j et k<>l}        if ((i=j) and (k<>l)) then begin

                              sumld^[i,k]:=1 + sumld^[i,k];
                              sumld^[i,l]:=1 + sumld^[i,l];
                              alefld^[i,lc2]:=2 + alefld^[i,lc2];
                              alefld^[k,lc1]:=1 + alefld^[k,lc1];
                              alefld^[l,lc1]:=1 + alefld^[l,lc1];
                              homold^[i,lc2]:=1 + homold^[i,lc2];
                              indld:= 1 + indld;

                              totsum^[i,k]:=1 +totsum^[i,k];
                              totsum^[i,l]:=1 +totsum^[i,l];
                              sumalf^[i,lc2]:=2 + sumalf^[i,lc2];
                              sumalf^[k,lc1]:=1 + sumalf^[k,lc1];
                              sumalf^[l,lc1]:=1 + sumalf^[l,lc1];
                              homosm^[i,lc2]:=1 + homosm^[i,lc2];
                              sumind:= 1 + sumind;
                            end;

{cas o— i<>j et k=l}        if ((i<>j) and (k=l)) then begin

                              sumld^[i,k]:=1 + sumld^[i,k];
                              sumld^[j,k]:=1 + sumld^[j,k];
                              alefld^[i,lc2]:=1 + alefld^[i,lc2];
                              alefld^[j,lc2]:=1 + alefld^[j,lc2];
                              alefld^[k,lc1]:=2 + alefld^[k,lc1];
                              homold^[k,lc1]:=1 + homold^[k,lc1];
                              indld:= 1 + indld;

                              totsum^[i,k]:=1 +totsum^[i,k];
                              totsum^[j,l]:=1 +totsum^[j,l];
                              sumalf^[i,lc2]:=1 + sumalf^[i,lc2];
                              sumalf^[j,lc2]:=1 + sumalf^[j,lc2];
                              sumalf^[k,lc1]:=2 + sumalf^[k,lc1];
                              homosm^[k,lc1]:=1 + homosm^[k,lc1];
                              sumind:= 1 + sumind;
                            end;

{cas o— i<>j et k<>l}       if ((i<>j) and (k<>l)) then begin

{on doit d‚clarer sumld en}   sumld^[i,k]:=0.5 + sumld^[i,k];
{single}                      sumld^[j,k]:=0.5 + sumld^[j,k];
                              sumld^[i,l]:=0.5 + sumld^[i,l];
                              sumld^[j,l]:=0.5 + sumld^[j,l];
                              alefld^[i,lc2]:=1 + alefld^[i,lc2];
                              alefld^[j,lc2]:=1 + alefld^[j,lc2];
                              alefld^[k,lc1]:=1 + alefld^[k,lc1];
                              alefld^[l,lc1]:=1 + alefld^[l,lc1];
                              indld:= 1 + indld;

{mˆme rem pour totsum que}    totsum^[i,k]:=0.5 +totsum^[i,k];
{pour sumld}                  totsum^[j,k]:=0.5 +totsum^[j,k];
                              totsum^[i,l]:=0.5 +totsum^[i,l];
                              totsum^[j,l]:=0.5 +totsum^[j,l];
                              sumalf^[i,lc2]:=1 + sumalf^[i,lc2];
                              sumalf^[j,lc2]:=1 + sumalf^[j,lc2];
                              sumalf^[k,lc1]:=1 + sumalf^[k,lc1];
                              sumalf^[l,lc1]:=1 + sumalf^[l,lc1];
                              sumind:= 1 + sumind;
                            end;
                          end; {condition ndr2[lc2]=l}
                        end;  {condition ndr2[lc2]<>0}
                      end; {compteur 2Šme allŠle 2Šme locus}
                    end; {condition ndr1[lc2]=k}
                  end; {condition ndr1[lc2]<>0}
                end; {compteur 1er allŠle 2Šme locus}
            end; {condition ndr2[lc1]=j}
          end;  {condition ndr2[lc1]<>0}
        end; {compteur 2Šme allŠle 1er locus}
      end; {condition ndr1[lc1]=i}
    end; {condition ndr1[lc1]<>0}
  end; {compteur 1er allŠle 1er locus}
 end; {boucle sur les individus}



        {boucle sur les sous-pop et les couples de
                    locus  dans le program principal}

end;

(**********************************************************************)
(**********   Calcul du desequilibre de liaison pop/pop   *************)
(**********************************************************************)


procedure ldcalc (allchk,l1,l2    :integer);

var       a1,a2,a3              :integer;
          ngfrc,ddl1,ddl2,d1,d2    :integer;
          sumbsq,sumde,sumz        :single;
          zstat,bsq,chi,zbar       :single;
          alefa,assoc,p,q,alefb    :single;
          bool                     :boolean;
          hwcora,hwcorb,denom      :single;

begin

       sumbsq:=0; sumde:=0; sumz:=0; ngfrc:=0; denom:=0;
       if (indld<>0) then begin
          bool:=false;

          for a1:=1 to nall[l1] do begin
              alefa:=alefld^[a1,l2]/(2*indld) ;
              if (alefa<>0) then begin

                 if (alefa=1) then begin

                  bool:=true;
                  for a2:=1 to nall[l2] do begin

                    assoc:=alefld^[a2,l1]/indld;
                    fone:=fone+sqr(assoc)*indld;
                    alefa:=alefld^[a2,l1]/(2*indld);
                    fseven:=fseven+sqr(alefa)*indld;

                    if (allchk=0) then begin
                      p:= alef[pop,l2,a2];
                      seven:=sqr(p)+seven;
                      fthree:=fthree+assoc*indld*alefa;
                    end;
                  end;

                 end {alefa=1}

                 else begin  {alefa < 1}

                  for a2:=1 to nall[l2] do begin

                    alefb:=alefld^[a2,l1]/(2*indld);

                    if (alefb<>0) then begin
 (*aucun des allŠles fixe*)
                      if (alefb < 1) then begin

                       assoc:=sumld^[a1,a2]/indld;
                       fone:=fone+sqr(assoc)*indld;
                       fthree:=fthree+assoc*indld*alefa*alefb;
                       fseven:=fseven+sqr(alefa)*sqr(alefb)*indld;

                       if (allchk=0) then begin
                         p:= alef[pop,l1,a1];
                         q:= alef[pop,l2,a2];
                         seven:=sqr(p)*sqr(q)+seven;
                       end;
                       hwcora:=(homold^[a1,l2]/indld)-sqr(alefa);
                       hwcorb:=(homold^[a2,l1]/indld)-sqr(alefb);
                       burrow^[a1,a2]:=(assoc-2*alefa*alefb)*(indld/(indld-1));
                       denom:=(alefa*(1-alefa)+hwcora)*
                                      (alefb*(1-alefb)+hwcorb);
                       if (denom <> 0) then begin
                            corrad^[a1,a2]:=burrow^[a1,a2]/sqrt(denom);
                            if (corrad^[a1,a2]>=1) then
                                      corrad^[a1,a2]:=0.99999;
                            if (corrad^[a1,a2]<(-1)) then
                               corrad^[a1,a2]:=(-0.99999);
                            zstat:=(0.5*Ln((1+corrad^[a1,a2])/
                                              (1-corrad^[a1,a2])));
                            zstat:=zstat*(indld-3);
                            sumz:=sumz+abs(zstat);
                            bsq:=sqr(abs(burrow^[a1,a2]));
                            chi:=indld*bsq/(alefa*alefb);
                            ngfrc:=ngfrc+1;
                            sumbsq:=sumbsq+bsq;
                            chisq:=chisq+chi;
                            sumde:=denom+sumde;
                       end; (* denom<>0*)
                      end              (*alefb <1*)
                      else begin       {alefb=1}
                       bool:=true;
                       for a3:=1 to nall[l1] do begin

                         assoc:=alefld^[a3,l2]/indld;
                         fone:=fone+sqr(assoc)*indld;
                         alefa:=alefld^[a3,l2]/(2*indld);
                         fseven:=fseven+sqr(alefa)*indld;

                         if (allchk=0) then begin
                            p:= alef[pop,l1,a3];
                            seven:=sqr(p)+seven;
                            fthree:=fthree+assoc*indld*alefa;
                         end;
                       end;
                      end;             {aleb=1}
                    end;           (*alefb<>0*)
                  end;          (*boucle sur a2*)
                 end;          (*alefa<1*)
              end;        {alefa<>0}
          end;         (*boucle sur a1*)



          if (bool=false)  then begin
            if (sumde<>0) then begin
              bursq:=sqrt(sumbsq/ngfrc);
              bursd:=sumbsq/sumde;
              if (indld-3)>0 then begin
              zbar:=sumz/((indld-3)*ngfrc);
              if exp(2*zbar)<>1 then begin
              z:=(exp(2*zbar)-1)/(exp(2*zbar)+1);

              {calcul des ddl pour les tests de chi2}

              ddl1:=0;
              ddl2:=0;

              for d1:=1 to nall[l1] do
              begin
                if alefld^[d1,l2]<>0 then begin
                   ddl1:=ddl1+1;
                end;
              end;
              for d2:=1 to nall[l2] do
              begin
                if alefld^[d2,l1]<>0 then begin
                   ddl2:=ddl2+1;
                end;
              end;

               negfre:=(ddl1-1)*(ddl2-1);
               df:=negfre; {def. du parametre df pour chiprb,outpop}
               prob:=(abs(1-chiprb(chisq,df)));
               if (prob<=0.0001) then prob:=0.0001;
               end; {condition exp(2*zbar<>1}
               end; {condition indld<>3}
            end;  {condition  sumde<>0}
          end    {bool=false}

          else begin
            fixchk:=1;
          end;    {condition 1 allŠle fix‚ … au moins un des 2 locus}

       end {condition indld[l1,l2]<>0}

       else begin        {un des locus est manquant}
          nodat:=1;
          npopct:=npopct-1
       end;
end; {ldcalc}



(******************************************************************)
(******R‚initialisations des matrices pour ldcalc pop/pop *********)
(******************************************************************)
procedure clear (a,b  :integer);

var  aa,bb      :integer;

begin

        indld:=0;
        fixchk:=0;
        nodat:=0;
        negfre:=0;
        bursq:=0;
        bursd:=0;
        chisq:=0;
        z:=0;
        prob:=0;
        for aa:=1 to nall[a] do
          for bb:=1 to nall[b] do begin
            sumld^[aa,bb]:=0;
            burrow^[aa,bb]:=0;
            corrad^[aa,bb]:=0;
            alefld^[aa,b]:=0;
            alefld^[bb,a]:=0;
            homold^[aa,b]:=0;
            homold^[bb,a]:=0;
            chisc[aa,bb]:=0;
            rcprob[aa,bb]:=0;
        end;

end;  {clear}

(******************************************************************)
(*transformation des matrices pour ldcalc dans toute la population*)
(******************************************************************)

procedure trans (c,d  :integer);

var cc,dd   :integer;

begin
        indld:=sumind;
        for cc:=1 to nall[c] do begin
          for dd:=1 to nall[d] do begin
            sumld^[cc,dd]:=totsum^[cc,dd];
            alefld^[cc,d]:=sumalf^[cc,d];
            alefld^[dd,c]:=sumalf^[dd,c];
            homold^[cc,d]:=homosm^[cc,d];
            homold^[dd,c]:=homosm^[dd,c];
          end;
        end;


end; {trans}

(******************************************************************)
(**********Calcul des fr‚quences all‚liques pop/pop****************)
(******************************************************************)

procedure alefp;

var   wgt                : ptstat;
      pp,ll,ic,ip     :integer;
      ipp,icc,ial        :integer;
      jpp,jcc,jal        :integer;

begin

    new(wgt);
    {initø de la matrice wgt}
    for pp:=1 to npop do
      for ll:=1 to nloc do begin
          wgt^[pp,ll]:=0;
    end;

    {calcul du nb total d'individus dans les populations}
    for ic:=1 to nloc do begin
      for ip:=1 to npop do begin
        totn[ic]:=totn[ic]+ind[ip,ic];
      end;
    end;

    {calcul des freq. all‚liques et idem pond‚r‚es}
    for ipp:=1 to npop do begin
      for icc:=1 to nloc do begin
        for ial:=1 to nall[icc] do begin
          if ind[ipp,icc]<>0 then begin
            alef[ipp,icc,ial]:=alef[ipp,icc,ial]/(2*ind[ipp,icc]);
          end;
        end;
      end;
    end;

    for jcc:=1 to nloc do
      for jal:=1 to nall[jcc] do
        for jpp:=1 to npop do begin
          if ind[jpp,jcc]<>0 then begin
            wgt^[jpp,jcc]:=ind[jpp,jcc]/totn[jcc];
            wta^[jcc,jal]:=wta^[jcc,jal]+(alef[jpp,jcc,jal]*wgt^[jpp,jcc])
          end;
        end;
     dispose(wgt);
end;

(******************************************************************)
(***********Edition du tableau des fr‚quences all‚liques***********)
(******************************************************************)

procedure freq(outfreq : string);

var   file4        :text;
      al,m   :integer;


begin

assign (file4,outfreq);
rewrite(file4);
    writeln (file4);
    writeln(file4,'        ALLELIC FREQUENCIES IN SUBPOPULATIONS        ');
    writeln (file4); writeln (file4);
    write (file4,'LOCUS  TOT');

   for m:=1 to 6 do begin
    pop:=(10*m-9);
    if pop <=npop then begin
     while (pop<=10*m) and (pop<=npop) do begin
         write (file4,'    ',pop,'  ');
         pop:=pop+1;
      end;
      writeln(file4);
      for loc:=1 to nloc do begin
         writeln (file4,aloc[loc]);
         write (file4,' (N)','  ',totn[loc]:4);
          pop:=(10*m-9);
          while (pop<=10*m) and (pop<=npop) do begin
           write (file4,'    ',ind[pop,loc]:3);
           pop:=pop+1;
          end;
          writeln(file4);
          for al:=1 to nall[loc] do begin
             write (file4,'  ',al,'  ',wta^[loc,al]:5:3);
             pop:=(10*m-9);
             while (pop<=10*m) and (pop<=npop) do begin
               write (file4,'  ',alef[pop,loc,al]:5:3);
               pop:=pop+1;
             end;  {boucle sur les populations}
             writeln(file4);
          end; {boucle sur les alleles}
        writeln(file4);
      end; {boucle sur les locus}
      writeln(file4);
    end; {condition sur pop}
   end; {boucle sur m}
   close (file4);

end; {freq}



(******************************************************************)
(*********    Decision si                                **********)
(***Ecriture detaillee des LD par couple d'alleles quand        ***)
(**** -le test global pour un couple de locus est significatif*****)
(**** -un couple d"alleles est fortement correle              *****)
(******************************************************************)

procedure checkprob (loci, locii  :integer; var choix :byte);

var  al1,al2,count     :integer;

begin

  count:=0;
  if (nodat=0) then
   begin
     for al1:=1 to nall[loci] do begin
       for al2:=1 to nall[locii] do begin
         if burrow^[al1,al2]<>0 then begin
            chisc[al1,al2]:=(abs(power(corrad^[al1,al2],2)))*indld;
            df:=1;
            rcprob[al1,al2]:=(abs(1-chiprb(chisc[al1,al2],df)));
            if (rcprob[al1,al2]<=0.0001) then rcprob[al1,al2]:=0.0001;
            if (rcprob[al1,al2]<=signi) then count:=count+1;
         end;
       end;
     end;
     if (count>=1) then choix:=1 else choix:=0;
   end;

 end;   {checkprob}



(******************************************************************)
(***Ecriture detaillee des LD par couple d'alleles quand        ***)
(**** -le test global pour un couple de locus est significatif*****)
(**** -un couple d"alleles est fortement correle              *****)
(******************************************************************)

procedure ldlwrt(loci,locii,allchk   :integer);

var   bb,al1,al2                         :integer;

begin


  writeln (filet);
  writeln (filet);
  for bb:=1 to 80 do
      write (filet,'*');
  writeln (filet);
  writeln (filet,'detailed output ',aloc[loci],' - ',aloc[locii]);
  writeln (filet);

if (allchk=0) and (nodat=0) then begin

writeln(filet,'population ',apop[pop]);
for bb:=1 to 80 do write (filet,'*');
writeln (filet);
writeln(filet);
writeln(filet,'Sample size= ',indld,' individuals');
writeln(filet);
write (filet,'   ALLELES     ','    T(IJ)','      D(IJ)  ','      R(IJ)  ');
writeln (filet,' CHI-SQUARE','   PROB.');
 for al1:=1 to nall[loci] do begin
    for al2:=1 to nall[locii] do begin
       if burrow^[al1,al2]<>0 then begin
         write (filet,'    ',al1,'  -  ',al2,'        ');
         write (filet,sumld^[al1,al2]:5:1,'    ');
         write (filet,burrow^[al1,al2]:8:5,'    ');
         write (filet,corrad^[al1,al2]:8:5,'    ',chisc[al1,al2]:7:2);
         writeln (filet,'    ',rcprob[al1,al2]:6:4);
       end;
    end; {boucle sur al2}
  end; {boucle sur al1}
writeln (filet);
writeln (filet,'variance= ',bursd:7:5);
writeln (filet,'common correlation= ', z:7:5);
writeln (filet,'chi-square(',negfre,')=',chisq:4:2);
writeln (filet,'prob.= ',prob:6:4);
end; {condition allchk=0 et nodat=0}

if (allchk=1) and (nodat=0) then begin

writeln (filet,'all subpopulations');
for bb:=1 to 80 do write (filet,'*');
writeln (filet);
writeln(filet);
writeln(filet,'Sample size= ',indld,' individuals');
writeln(filet);
write (filet,'   ALLELES     ','    T(IJ)','      D(IJ)  ','      R(IJ)  ');
writeln (filet,' CHI-SQUARE','   PROB.');
 for al1:=1 to nall[loci] do begin
    for al2:=1 to nall[locii] do begin
       if burrow^[al1,al2]<>0 then begin
         write (filet,'    ',al1,'  -  ',al2,'        ');
         write (filet,sumld^[al1,al2]:5:1,'    ');
         write (filet,burrow^[al1,al2]:8:5,'    ');
         write (filet,corrad^[al1,al2]:8:5,'    ',chisc[al1,al2]:7:2);
         writeln (filet,'    ',rcprob[al1,al2]:6:4);
       end;
    end; {boucle sur al2}
  end; {boucle sur al1}
writeln (filet);
writeln (filet,'variance= ',bursd:7:5);
writeln (filet,'common correlation= ', z:7:5);
writeln (filet,'chi-square(',negfre,')=',chisq:4:2);
writeln (filet,'prob.= ',prob:6:4);
end; {condition allchk=1 et nodat=0}


end; {ldlwrt}




(******************************************************************)
(***************Ecriture condensee des LD**************************)
(******************************************************************)

procedure outpop (allchk  :integer;signi  :single;loci,locii  :integer);

var       bb,sig          :integer;
          choose          :byte;

begin


       if (fixchk>=1) then begin
          if (allchk=0) then write (file2,apop[pop],'             ')
             else write (file2,'ALL SUBPOP       ');
          writeln (file2,': AN ALLELE AT ONE OR BOTH LOCI IS FIXED')
       end
       else if (nodat=1) then begin
              if (allchk=0) then write (file2,apop[pop],'             ')
                  else write (file2,'ALL SUBPOP       ');
              writeln (file2,'NO INFORMATION WAS COLLECTED UPON ONE OF THE LOCUS')
            end
       else begin
           if (allchk=0) then write (file2,apop[pop],'             ')
                            else write (file2,'ALL SUBPOP       ');
           write (file2,'    ',indld:3,'        ');
           write (file2,'  ',z:7:5,'      ');
           write (file2,' ',chisq:4:2,'     ');
           writeln (file2,negfre,'    ',prob:6:4)
       end;
    sig:=0;
       if (fixchk=0) then
          if (nodat=0) then begin
             checkprob (loci,locii,choose);
             if ((prob<=signi) or (choose=1)) then begin
                sig:=1;
                ldlwrt(loci,locii,allchk)
             end;
          end;
end;  {outpop}



(******************************************************************)
(*********Ecriture des composantes de la variance du LD************)
(******************************************************************)

procedure ohta (loc1,loc2  :integer);

var       twelve,ftwo,fsix,fthirt,thirt   :single;
          dit,dstp,disp,dis,dst           :single;
          xbar,ybar,gbar                  :single;
          fthre,fon,fsevn,sevn,diff       :single;
          aa,al1,al2                      :integer;
          file4                           :text;


begin
  
       if (npopct > 1) then begin
          twelve:=0; ftwo:=0; fsix:=0;
          fthirt:=0; thirt:=0;
          for al1:=1 to nall[loc1] do
           for al2:=1 to nall[loc2] do begin
              xbar:=sumalf^[al1,loc2]/(2*sumind);
              ybar:=sumalf^[al2,loc1]/(2*sumind);
              gbar:=totsum^[al1,al2]/sumind;
              ftwo:=ftwo+sqr(gbar);
              fsix:=fsix+2*xbar*ybar*gbar;
              fthirt:=fthirt+4*sqr(xbar)*sqr(ybar);
              thirt:=thirt+4*sqr(wta^[loc1,al1])*sqr(wta^[loc2,al2]);
              for pop:=1 to npop do begin
               if ((alef[pop,loc1,al1]<=1) and (alef[pop,loc2,al2]<=1))
                  then begin
                         diff:=alef[pop,loc1,al1]*alef[pop,loc2,al2]
                                   *wta^[loc1,al1]*wta^[loc2,al2];
                         twelve:=twelve+diff;
                  end;  {alleles inf a 1}
              end;      {boucle populations}
          end;      {boucle alleles}

          fthre:=(2*fthree)/sumind;
          fon:=fone/sumind;
          fsevn:=4*fseven/sumind;
          sevn:=4*seven/npopct;
          twelve:=4*twelve/npopct;
          dit:=fon+fthirt-2*fsix;
          dstp:=ftwo+fthirt-2*fsix;
          disp:=fon-ftwo;
          dis:=fon+fsevn-2*fthre;
          dst:=sevn+thirt-2*twelve;
          dst:=dst/4;

          if ((dit<=0) and (dit>= -0.000001)) then dit:=0;
          if ((dis<=0) and (dis>= -0.000001)) then dis:=0;
          if ((disp<=0) and (disp>= -0.000001)) then disp:=0;
          if ((dst<=0) and (dst>= -0.000001)) then dst:=0;
          if ((dstp<=0) and (dstp>= -0.000001)) then dstp:=0;

          write (file3,'  ',aloc[loc1],' - ',aloc[loc2],'    ');
          write (file3,'  ',dis:7:5,'  ',disp:7:5,'      ');
          write (file3,'  ',dst:7:5,'  ',dstp:7:5,'      ');
          writeln (file3,'   ',dit:7:5);

       end;    {npopct superieur a 1}


end;  {ohta}



(*****************************************************************)
(**********  En tetes des fichiers de sortie   *******************)
(*****************************************************************)


(**********Ecriture condensee des LD ***********)


procedure condldwrt (lc1,lc2  :integer);

var   a      :integer;

begin

 writeln (file2); writeln (file2);
 for a:=1 to 80 do write (file2,'*');
 writeln (file2);
 writeln (file2,'ANALYSIS OF LINKAGE DESEQUILIBRIUM BETWEEN ',aloc[lc1],' - ',aloc[lc2]);
 writeln (file2);
 for a:=1 to 80 do
      write (file2,'*');
  writeln (file2);
  writeln (file2);
  writeln (file2,'                 NUMBER OF      COMMON         CHI-');
  write (file2,'  POPULATION      COMPARISONS    CORRELATION    SQUARE');
  writeln (file2,'    D.F.    PROB.');
  for a:=1 to 80 do
      write (file2,'-');
  writeln (file2);

 end;   {condldwrt}



(************** Ecriture des composantes de la variance du deseq ***)


procedure ohtawrt;

var  aa          :integer;

begin

 for aa:=1 to 80 do
      write (file3,'*');
  writeln (file3);
  writeln (file3,'VARIANCE COMPONENTS OF LINKAGE DESEQUILIBRIUM');
  for aa:=1 to 80 do
      write (file3,'-');
  writeln (file3);
  write (file3,'                   WITHIN SUBPOPULATION    ');
  write (file3,'BETWEEN SUBPOPULATION    TOTAL POPULATION');
  writeln (file3);
  writeln (file3,'LOCI COMPARATED         COMPONENTS              COMPONENT');
  write (file3,'                   --------------------    ');
  writeln (file3,'--------------------    ----------------');
  write (file3,'                     D(IS)      D"(IS)       D(ST)      D"(ST)');
  writeln (file3,'          D(IT)');
  for aa:=1 to 80 do
      write (file3,'-');
  writeln (file3);

end;  {ohtawrt}


(******************************************************************)
(*********************** affichage de messages a l"ecran***********)
(******************************************************************)

procedure screen1;

var  aa          :integer;
begin
  clrscr;
  for aa:=1 to 12 do
     writeln;
  writeln ('               COMPUTING ALLELIC FREQUENCIES ...');
end;


procedure screen2 (l1,l2  :integer);

var  aa           :integer;

begin
  clrscr;
  for aa:=1 to 12 do writeln;
  writeln ('COMPUTING LINKAGE DESEQUILIBRIA BETWEEN ',aloc[l1],'  -  ',aloc[l2]);
end;



(******************************************************************)
(***********************Programme principal************************)
(******************************************************************)


begin

      Clrscr;
      writeln ('Please answer Y or N to the following questions:');
      writeln ('Do you want the LD computation subpop/subpop?');
      readln (response);
      if ((response='Y') or (response='y')) then ldopt:=1
                                            else ldopt:=0;
      writeln ('Do you want the LD computation in all subpopulations?');
      readln (response);
      if ((response='Y') or (response='y')) then ldtotopt:=1
                                            else ldtotopt:=0;

      if ((ldtotopt=1) or  (ldopt=1))
         then begin
           writeln ('What Significance level do you wish?(0.01,0.05,or 0.1 for ex):');
           readln (signi)
         end;
      writeln ('Do you want the Ohta"s LD components variance analysis');
      writeln ('to be performed?');readln (response);
      if ((response='Y') or (response='y')) then ohtaopt:=1
                                            else ohtaopt:=0;

      writeln ('Do you want an intrasubpopulation allele frequencies table?');
      readln (response);
      if ((response='Y') or (response='y')) then freqopt:=1
                                            else freqopt:=0;

      writeln ('Input data file-Populations names,loci names,');
      write ('Genotypes coding-:'); readln (Input);
      if ((ldopt=1)  or (ldtotopt=1))
      then begin
        write ('LD output file:'); readln (outdes)
      end;
      if (freqopt=1) then begin
      write ('Allelic frequencies output file:'); readln (outfreq)
      end;
      if (ohtaopt=1) then begin
      write ('Ohta LD variance components file:'); readln (outohta)
      end ;

          new(alefld);new(homold);new(sumalf);new(homosm);
          new(sumld);new(totsum);new(burrow);new(corrad);
          new(wta);


      init;

      assign (file1,input);
      reset (file1);
     {lecture du nombre de pop et du nombre de locus}
        read (file1,npop,nloc);
        readln (file1);
     {lecture des noms de populations et du nb d'ind/population}
        for pop:=1 to npop do
           begin
             readln (file1,apop[pop]);
             readln (file1,nmax[pop]);
           end;
      {lecture des noms des locus et du nb d'allŠles/locus}
        for loc:=1 to nloc do
           begin
             read (file1,aloc[loc],nall[loc]);readln (file1);
           end;
      {lecture des g‚notypes par pop,locus et couple de locus}
        for pop:=1 to npop do begin
            lecture;
        end;

        screen1;
        alefp;          {calcul des frequences alleliques}

        if (freqopt=1) then begin
          freq(outfreq);
        end;

        close (file1);

      if ((ldopt=1) or (ldtotopt=1)) then begin
         assign (file2,outdes);
         rewrite (file2);
      end;                          {ouverture des fichiers}
      if (ohtaopt=1) then begin
         assign (file3,outohta);
         rewrite (file3);
         ohtawrt;
      end;
      assign (filet,'tempo');


      for loc:=1 to (nloc-1) do
       for loco:=(loc+1) to nloc do
         begin
                                       {boucle sur les couples de locus}

           screen2 (loc,loco);
           initcouple (loc,loco);
           if ((ldopt=1) or (ldtotopt=1)) then begin
              rewrite (filet);
              condldwrt(loc,loco);
           end;

           reset (file1);
           for cc:=1 to (1+2*npop+nloc) do
                                 readln (file1);

           for pop:=1 to npop do
               begin
                 allchk:=0;
                 remp(loc,loco);
                 if ((ldopt=1) or (ohtaopt=1)) then begin
                    ldcalc(allchk,loc,loco);
                    if (ldopt=1) then outpop (allchk,signi,loc,loco);
                 end;

                 clear (loc,loco);       {reinitialisation des matrices }
                                         {pour le test du chi2}
               end;  {boucle sur les pop}

            close (file1);

            if (ohtaopt=1) then begin
                allchk:=0;
                ohta (loc,loco);
            end;                         {ohta}

            if (ldtotopt=1) then begin
                allchk:=1;
                trans (loc,loco);  {transformation des matrices indld en sumind et alefld en sumalf}
                ldcalc (allchk,loc,loco);
                outpop(allchk,signi,loc,loco);
            end;

            if ((ldopt=1) or (ldtotopt=1)) then begin
                close (filet);
                reset (filet);
                while not eof(filet) do begin
                 read (filet,response);
                 write (file2,response);
                end;
               close (filet);
            end;                               {population totale}

         end;                             {boucle sur les locus}


      if ((ldopt=1) or (ldtotopt=1)) then close (file2);
      if (ohtaopt=1) then close (file3);


  dispose(alefld);dispose(homold);dispose(sumalf);dispose(homosm);
  dispose(sumld);dispose(totsum);dispose(burrow);dispose(corrad);
  dispose(wta);

end.
