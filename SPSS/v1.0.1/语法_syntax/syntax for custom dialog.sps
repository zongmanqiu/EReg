* Encoding: UTF-8.
PRESERVE.
set printback=off.
set LENGTH=256.
set mxloops=1000000000.
/** standardized coefficient for polynomial model **/.
define @polystdb(@polyb=!charend('/')/@polysd=!charend('/')/@polym=!charend('/'))
compute @polytb=!@polyb.
do if (!@polym=6).
    compute @polytb(4,1)=@polytb(4,1)/!@polysd(4,1)*!@polysd(2,1)*!@polysd(3,1).
    else if (!@polym=7).
    compute @polytb(3,1)=@polytb(3,1)/!@polysd(3,1)*!@polysd(2,1)*!@polysd(2,1).
    compute @polytb(5,1)=@polytb(5,1)/!@polysd(5,1)*!@polysd(2,1)*!@polysd(4,1).
    else if (!@polym=8).
    compute @polytb(3,1)=@polytb(3,1)/!@polysd(3,1)*!@polysd(2,1)*!@polysd(2,1).
    compute @polytb(5,1)=@polytb(5,1)/!@polysd(5,1)*!@polysd(2,1)*!@polysd(4,1).
    compute @polytb(6,1)=@polytb(6,1)/!@polysd(6,1)*!@polysd(2,1)*!@polysd(2,1)*!@polysd(4,1).
    else if (!@polym=9).
    compute @polytb(4,1)=@polytb(4,1)/!@polysd(4,1)*!@polysd(3,1)*!@polysd(3,1).
    compute @polytb(5,1)=@polytb(5,1)/!@polysd(5,1)*!@polysd(2,1)*!@polysd(3,1). 
    else if (!@polym=10).
    compute @polytb(4,1)=@polytb(4,1)/!@polysd(4,1)*!@polysd(3,1)*!@polysd(3,1).
    compute @polytb(5,1)=@polytb(5,1)/!@polysd(5,1)*!@polysd(2,1)*!@polysd(3,1). 
    compute @polytb(6,1)=@polytb(6,1)/!@polysd(6,1)*!@polysd(2,1)*!@polysd(3,1)*!@polysd(3,1).
    else if (!@polym=11).
    compute @polytb(3,1)=@polytb(3,1)/!@polysd(3,1)*!@polysd(2,1)*!@polysd(2,1).
end if.
compute !@polyb=@polytb.
release @polytb.
!enddefine.
/** pai, product of a sequence**/.
define @Cpai (@Craw=!charend('/'))
    compute @Cr=nrow(!@Craw).
    compute @Cc=ncol(!@Craw).
    compute @Cpaiv=1.
    loop #cpai1=1 to @Cr.
        loop #cpai2= 1 to @Cc.
            compute @Cpaiv=@Cpaiv*!@Craw(#cpai1,#cpai2).
        end loop.
    end loop.
release @Cr,@Cc.
!enddefine.
/** Ridge Regression version 3 **/. 
define @RR3(@RRyx=!charend('/')/@RRk=!default(0) !charend('/')/@RRext=!default(0) !charend('/')
/@RRkest=!default(0) !charend('/'))
compute @RRyx=!@RRyx.
compute @RRk=!@RRk.
compute @RRnc=nrow(@RRyx).
compute @RRnx=ncol(@RRyx)-1.
compute @RRnk=ncol(@RRk).
/** variance-covariance matrix of variables **/.
compute @RRvcv=(sscp(@RRyx)-(sscp(csum(@RRyx))/@RRnc))/(@RRnc-1).
compute @RRsd=sqrt(diag(@RRvcv)).  /** std.dev **/.
compute @RRyxr=inv(mdiag(@RRsd))*@RRvcv*inv(mdiag(@RRsd))./** correlation matrix **/.
compute @RRyxm=csum(@RRyx)/@RRnc.  /** mean **/.
compute @RRyx=@RRyx-make(@RRnc,1,1)*@RRyxm. /** mean center **/.
compute @RRx=@RRyx(:,2:(@RRnx+1)).
compute @RRy=@RRyx(:,1).
compute @RRxsc=csum(@RRx&**2)&**0.5. /** scaling factor **/.
/** scaling x in correlation form **/.
compute @RRx=@RRx&/(make(@RRnc,1,1)*@RRxsc). 
call svd(@RRx,@RRsvdu,@RRsvdd,@RRsvdv).
compute @RRsvdu=@RRsvdu(:,1:@RRnx).
compute @RRsvdd=diag(@RRsvdd(1:@RRnx,1:@RRnx)).
compute @RRrhs=t(@RRsvdu)*@RRy.
compute @RRdiv=(@RRsvdd*make(1,@RRnk,1))&**2+make(@RRnx,1,1)*@RRk.
compute @RRb=diag(@RRsvdd*t(@RRrhs))*make(1,@RRnk,1)&/@RRdiv.
loop #rrnk1=1 to @RRnk.
    compute @RRb(:,#rrnk1)=@RRsvdv*@RRb(:,#rrnk1).
end loop.
compute @RRbr=@RRb. /** ridge coefficients **/.
compute @RRb=@RRb&/(t(@RRxsc)*make(1,@RRnk,1)). /** raw coefficients **/.
/** stdyx.coefficients **/.
compute @RRbstd=@RRb/@RRsd(1,1)&*(@RRsd(2:(@RRnx+1),1)*make(1,@RRnk,1)).
/** intercepts **/.
compute @RRb={@RRyxm(1,1)-@RRyxm(1,2:(@RRnx+1))*@RRb;@RRb}.
do if (!@RRext=1).
    compute @RRrfit=@RRx*@RRbr.    /** ridge fit **/.
    compute @RRres=@RRy*make(1,@RRnk,1)-@RRrfit. /** ridge residuals **/.
    compute @RRsst=sscp(@RRy)*make(1,@rrnk,1). /** sum of square total **/.
    compute @RRssr=csum(@RRrfit&**2).  /** sum of square regression **/.
    compute @RRsse=csum(@RRres&**2). /** sum of square error **/.
    compute @RRrsq=@RRssr&/@RRsst.    /** r square **/.
    /** adjusted r square **/.
    compute @RRrsqa=1-(1-@RRrsq)*((@RRnc-1)/(@RRnc-@RRnx)).  
    compute @RRxr=@RRyxr(2:(@RRnx+1),2:(@RRnx+1)).    /** predictors correlation matrix **/.
    compute @RRident=ident(@RRnx,@RRnx).
    compute @RRres2s=csum(@RRres&**2).
    compute @RRbrse=make(@RRnx,@RRnk,0). /** standard error of ridge coefficients**/.
    compute @RRb0se=make(1,@RRnk,0).         /** standard error of intercept **/.
    compute @RRb0=make(1,@RRnk,0).             /** intercenpt **.
    compute @RRdft=make(1,@RRnk,0).             /** degreed of freedom of coefficients **/.
    compute @RRdf2=make(1,@RRnk,0).            /** df2 of model F test, **/.
    compute @RRFv=make(1,@RRnk,0).             /** F value of model **.
    compute @RRFp=make(1,@RRnk,0).            /** p value of F **/.
    compute @RRaic=make(1,@RRnk,0).            /** AIC **/.
    compute @RRbic=make(1,@RRnk,0).            /** BIC **/.
    compute @RRmse=make(1,@RRnk,0).          /** MSE **/.
    compute @RRcv=make(1,@RRnk,0).             /** Cross Validation, CV **/.
    compute @RRvif=make(@RRnx,@RRnk,0).    /** VIF **/.
    @reg3 @rgy=@RRy/@rgx=@RRx/@rgci=95 .    /** OLS regression **/.
    compute rgb=rgb(2:(@RRnx+1),1).     /** OLS coefficients **/.
    call eigen(@RRxr,@RRvec,@RRval).  /** eigenvalue and eigenvector **/.
    loop #nk1=1 to @RRnk.
        compute @RRtmp0=inv(@RRxr+@RRident*@RRk(1,#nk1)).
        compute @RRvif(:,#nk1)=diag(@RRtmp0*@RRxr*@RRtmp0). /** vif **/.
        compute @RRtmp1=@RRtmp0*t(@RRx).   /** z **/.
        compute @RRtmp2=@RRx*@RRtmp1.      /** hatr **/.
        compute @RRcv(1,#nk1)=csum((@RRres(:,#nk1)/
            (1-(1/@RRnc)-diag(@RRtmp2)))&**2)/@RRnc.
        compute @RRdft(1,#nk1)=@RRnc-csum(diag(@RRtmp2)).   /** b df **/.
        compute @RRtmp3=@RRnc-csum(diag(2*@RRtmp2-sscp(t(@RRtmp2)))).   /** redf **/.
        compute @RRdf2(1,#nk1)=@RRtmp3.
        compute @RRtmp5=@RRres2s(1,#nk1)/@RRtmp3.   /** rsigma2 **/.
        /** variance-covariance matrix of coefficients **/.
        compute @RRtmp6=sscp(t(@RRtmp1))*@RRtmp5.   /** vcov **/.
        compute @RRbrse(:,#nk1)=sqrt(diag(@RRtmp6)). /** ridge standard errors **/.
        compute @RRb0se(:,#nk1)=sqrt(@RRsd(1,1)**2/
            @RRnc+@RRyxm(1,2:(@RRnx+1))&**2*(@RRbrse(:,#nk1)&**2)). /** ridge b0se **/.
        /** ridge b0 **/.
        compute @RRb0(1,#nk1)=@RRyxm(1,1)-@RRyxm(1,2:(@RRnx+1))*@RRbr(:,#nk1).  
        compute @RRFv(1,#nk1)=1/@RRnx*(t(@RRbr(:,#nk1))*
            inv(@RRtmp6)*@RRbr(:,#nk1) ).  /** F value of model **/.  
        compute @RRtmp7=csum((-@RRk(1,#nk1)*@RRtmp0*rgb)&**2).   /** bias **/.
        compute @RRtmp8=csum(@RRval&/(@RRval+@RRk(1,#nk1))&**2).  /** rvarcal **/.
        compute @RRmse(1,#nk1)=@RRtmp8*@RRtmp5+@RRtmp7.
    end loop.
    compute @RRdf1=@RRnc-@RRdft. /** df1 of model F test, **/.
    compute @RRdft=@RRdft-1.
    compute @RRbr={@RRb0;@RRbr}.
    compute @RRbrse={@RRb0se;@RRbrse}.
    compute @RRbt=@RRbr&/@RRbrse.    /** coefficients t value **/.
    compute @RRbp=@RRbrse.    /** coefficients p value **/.
    loop #nk2=1 to @RRnk. 
        compute @RRaic(1,#nk2)=@RRnc*ln(@RRsse(1,#nk2)/@RRnc)+2*@RRdf1(1,#nk2).
        compute @RRbic(1,#nk2)=@RRnc*ln(@RRsse(1,#nk2))+@RRdf1(1,#nk2)*ln(@RRnc).
        compute @RRFp(1,#nk2)=1-fcdf(@RRFv(1,#nk2),@RRdf1(1,#nk2),@RRdf2(1,#nk2)).
        compute @RRbp(:,#nk2)=2*(1-tcdf(abs(@RRbt(:,#nk2)),@RRdft(1,#nk2))).
    end loop.
    compute @RRbse=@RRb&/@RRbt.
    compute @RRmfit={@RRrsq;@RRrsqa;@RRmse;@RRaic;@RRbic}.
    compute @RRmtest={@RRFv;@RRdf1;@RRdf2;@RRFp;@RRdft}.
    do if (!@RRkest=1).
        compute @RRpred=@RRx*rgb.
        compute @RRsgma=csum((@RRy-@RRpred)&**2)/(@RRnc-@RRnx).
        compute @RRxt=@RRx*@RRvec.
        compute @RRahat=inv(mdiag(@RRval))*t(@RRxt)*@RRy.
        /** HK 1970 **/.
        compute @RRks=@RRsgma/cmax(@RRahat&**2)..
        /** HKB 1975 **/.
        compute @RRks={@RRks;@RRnx*@RRsgma/csum(rgb&**2)}.
        /** LW 1976 **/.
        compute @RRks={@RRks;@RRnx*@RRsgma/csum(@RRval&*@RRahat&**2)}. 
        /** HSL 1976 **/.
        compute @RRks={@RRks;@RRsgma*csum(@RRval&*@RRahat)**2/
            (csum(@RRval&*@RRahat&**2))**2}.
        /** AM 2003 **/.
        compute @RRks={@RRks;1/@RRnx*csum(@RRsgma&/(rgb&**2))}.
        /** GM 2003 **/.
        @Cpai @Craw={@RRahat&**2}.
        compute @RRks={@RRks;@RRsgma/(@Cpaiv&**(1/@RRnx))}.
        /** MED 2003 **/.
        compute @RRtmp1={@RRsgma&/(@RRahat&**2)}.
        @sortv sortv=@RRtmp1/sorttp=1.
        @prcn prcny=@RRtmp1/prcnp=50.
        compute @RRks={@RRks;percnv}.
        /** KS 2005 **/.
        compute @RRks={@RRks;cmax(@RRval)*@RRsgma/
            ((@RRnc-@RRnx)*@RRsgma+cmax(@RRval)*cmax(@RRsgma&/(@RRahat&**2)))}.
        compute @RRtmp1=@RRval*@RRsgma&/((@RRnc-@RRnx)*
            @RRsgma+@RRval&*@RRahat).
        /** KSarith 2005 **/.
        compute @RRks={@RRks;csum(@RRtmp1)/@RRnx}.
        /** KSmax 2005 **/.
        compute @RRks={@RRks;cmax(@RRtmp1)}.
        /** KSmd 2005 **/.
        @sortv sortv=@RRtmp1/sorttp=1.
        @prcn prcny=@RRtmp1/prcnp=50.
        compute @RRks={@RRks;percnv}.
        /** M1 2009 **/.
        compute @RRtmp1=@RRval*@RRsgma&/
            (@RRval&*(@RRahat&**2)+(@RRnc-@RRnx)*@RRsgma).
        @Cpai @Craw=@RRtmp1.
        compute @RRks={@RRks;@Cpaiv&**(1/@RRnx)}.
        compute @RRtmp1=sqrt(@RRsgma/@RRahat&**2).
        /** M2 2009 **/.
        compute @RRks={@RRks; cmax(1/@RRtmp1)}.
        /** M3 2009 **/.
        compute @RRks={@RRks; cmax(@RRtmp1)}.
        /** M4 2009 **/.
        compute @RRtmp2=1/@RRtmp1.
        @Cpai @Craw=@RRtmp2.
        compute @RRks={@RRks;@Cpaiv**(1/@RRnx)}.
        /** M5 2009 **/.
        compute @RRks={@RRks;(1/@Cpaiv)**(1/@RRnx)}.
        /** M6 2009 **/.
        @sortv sortv=@RRtmp2/sorttp=1.
        @prcn prcny=@RRtmp2/prcnp=50.
        compute @RRks={@RRks;percnv}.
        /** M7 2009 **/.
        @sortv sortv=@RRtmp1/sorttp=1.
        @prcn prcny=@RRtmp1/prcnp=50.
        compute @RRks={@RRks;percnv}.
        compute @RRtmp1=sqrt(cmax(@RRval)*@RRsgma/
            ((@RRnc-@RRnx)*@RRsgma+@RRahat&**2*cmax(@RRval))).
        /** M8 2012 **/.
        compute @RRks={@RRks;cmax(1/@RRtmp1)}.
        /** M9 2012 **/. 
        compute @RRks={@RRks;cmax(@RRtmp1)}.
        /** M10 2012 **/. 
        @Cpai @Craw=@RRtmp1.
        compute @RRks={@RRks;(1/@Cpaiv)**(1/@RRnx)}.
        /** M11 2012 **/. 
        compute @RRks={@RRks;@Cpaiv**(1/@RRnx)}.
        /** M12 2012 **/. 
        compute @RRtmp2=1/@RRtmp1.
        @sortv sortv=@RRtmp2/sorttp=1.
        @prcn prcny=@RRtmp2/prcnp=50.
        compute @RRks={@RRks;percnv}.
        /** D 2010 **/.
        compute @RRtmp1={0,@RRks(2,1)-1/@RRnc*cmax(rmax(@RRvif))}.
        compute @RRks={@RRks;rmax(@RRtmp1)}.
        compute @RRtmp1=@RRahat&**2.
        /** 4AD 2014 **/.
        compute @RRtmp2={0,2*@RRnx/cmax(@RRval)*csum(@RRsgma&/@RRtmp1)}.
        compute @RRks={@RRks;rmin(@RRtmp2)}.
        /** CV and GCV **/.
        compute @RRgcv=make(1,@RRnk,0).    /** Generalized CV **/.
        loop #rrnk3=1 to @RRnk.
            compute @RRgcv(1,#rrnk3)=csum((@RRy-@RRx*@RRbr(2:(@RRnx+1),#rrnk3))
                &**2)/(@RRnc-csum(@RRsvdd&**2&/@RRdiv(:,#rrnk3)))&**2.
        end loop.
        compute @RRtmp3=grade(@RRcv).        
        compute @RRtmp4=grade(@RRgcv). 
        loop #rrnk4=1 to @RRnk.
            do if (@RRtmp3(1,#rrnk4)=1).
                compute @RRtmp5=@RRk(1,#rrnk4).
            end if.
            do if (@RRtmp4(1,#rrnk4)=1).
                compute @RRtmp6=@RRk(1,#rrnk4).
            end if.
        end loop.
        compute @RRks={@RRks;@RRtmp5;@RRtmp6}.
        /** rank **/.
        compute @RRtmp2=grade(-@RRks).
        compute @RRks={@RRks,@RRtmp2}.
        compute @RRksrnm={'K.HK','K.HKB','K.LW','K.HSL','K.AM','K.GM','K.MED',
            'K.KS','KSarith','K.KSmax','K.KSmd','K.M1','K.M2','K.M3','K.M4','K.M5','K.M6','K.M7',
            'K.M8','K.M9','K.M10','K.M11','K.M12','K.D','K.4AD','CV','GCV','       .'}.
        compute @RRkscnm={"K value","Rank"}.
    end if.
    release @RRtmp0,@RRtmp1,@RRtmp2,@RRtmp3,@RRtmp5,
        @RRtmp6,@RRtmp7,@RRvec,@RRval.
end if.
    release @RRyx,@RRx,@RRy,@RRsvdu,@RRsvdd,@RRsvdv,@RRdiv.
!enddefine. 
/** delete the number which is out of range **/.
define @num.out(@no=!charend('/')/@nmm=!charend('/'))
    compute @no=!@no.
    compute @nolp=nrow(@no).
    compute @noz=0.
    loop #no = 1 to @nolp.
        do if (@no(#no,1)>= cmin(!@nmm) and @no(#no,1)<= cmax(!@nmm)).
            compute @noz=@noz+1.
            compute !@no(@noz,1)=@no(#no,1).
        end if.
    end loop.
    compute !@no=!@no(1:@noz,1).
    release @no,@nolp,@noz.
!enddefine.
/** solve quadratic equation a*x^2+b*x+c=0 **/.
define @eq2.slv(@eq2.ipt=!charend('/'))
    compute @eq2.a=!@eq2.ipt(1,1).
    compute @eq2.b=!@eq2.ipt(1,2).
    compute @eq2.c=!@eq2.ipt(1,3).
    compute @eq2.dt=@eq2.b**2-4*@eq2.a*@eq2.c.
    do if (@eq2.dt>=0).
        do if (@eq2.dt=0).
            compute @eq.rt={-@eq2.b}.
            else.
            compute @eq.rt={-@eq2.b-sqrt(@eq2.dt);-@eq2.b+sqrt(@eq2.dt)}.
        end if.
        compute @eq.rt=@eq.rt/@eq2.a*0.5.
        compute @eq.rtc=1.
        else.
        compute @eq.rtc=0.
    end if.
    release @eq2.a,@eq2.b,@eq2.c,@eq2.dt.
!enddefine.
/** solve cubic equation a*x^3+b*x^2+c*x+d=0 **/.
define @eq3.slv (@eq3.ipt=!charend('/'))
    compute @eq3.a=!@eq3.ipt(1,1).
    compute @eq3.b=!@eq3.ipt(1,2).
    compute @eq3.c=!@eq3.ipt(1,3).
    compute @eq3.d=!@eq3.ipt(1,4).
    compute @eq3.a2=@eq3.b**2-3*@eq3.a*@eq3.c.
    compute @eq3.b2=@eq3.b*@eq3.c-9*@eq3.a*@eq3.d.
    compute @eq3.c2=@eq3.c**2-3*@eq3.b*@eq3.d.
    compute @eq3.dt=@eq3.b2**2-4*@eq3.a2*@eq3.c2.
    do if (@eq3.a2=0 and @eq3.b2=0).
        compute @eq.rt=-@eq3.b/(3*@eq3.a).
        else.
        do if (@eq3.dt>0).
            compute @eq3.t1=@eq3.a2*@eq3.b+1.5*@eq3.a*(-@eq3.b2-sqrt(@eq3.dt)).
            compute @eq3.t2=@eq3.a2*@eq3.b+1.5*@eq3.a*(-@eq3.b2+sqrt(@eq3.dt)).
            do if (@eq3.t1>=0).
                compute @eq3.t1=@eq3.t1**(1/3).
                else.
                compute @eq3.t1=-((-@eq3.t1)**(1/3)).
            end if.
            do if (@eq3.t2>=0).
                compute @eq3.t2=@eq3.t2**(1/3).
                else.
                compute @eq3.t2=-((-@eq3.t2)**(1/3)).
            end if.
            compute @eq.rt=(-@eq3.b-@eq3.t1-@eq3.t2)/(3*@eq3.a).
            release @eq3.t1,@eq3.t2.
            else if (@eq3.dt=0).
            compute @eq3.t1=@eq3.b2/@eq3.a2.
            compute @eq.rt={-@eq3.b/@eq3.a+@eq3.t1;-@eq3.t1/2}.
            release @eq3.t1.
            else if (@eq3.dt<0).
            compute @eq3.t1=(2*@eq3.a2*@eq3.b-3*@eq3.a*@eq3.b2)/(2*sqrt(@eq3.a2**3)).
            compute @eq3.t2=2*artan(1)-arsin(@eq3.t1).
            compute @eq.rt=(-@eq3.b-2*sqrt(@eq3.a2)*cos(@eq3.t2/3))/(3*@eq3.a).
            compute @eq.rt={@eq.rt;(-@eq3.b+sqrt(@eq3.a2)*(cos(@eq3.t2/3)
                -sqrt(3)*sin(@eq3.t2/3)))/(3*@eq3.a)}.
            compute @eq.rt={@eq.rt;(-@eq3.b+sqrt(@eq3.a2)*(cos(@eq3.t2/3)+sqrt(3)
                *sin(@eq3.t2/3)))/(3*@eq3.a)}.
            release @eq3.t1,@eq3.t2.
        end if.
    end if.
    release @eq3.a,@eq3.b,@eq3.c,@eq3.d,@eq3.a2,@eq3.b2,@eq3.c2,@eq3.dt.
!enddefine.
/** solve quartic equation a*x^4+b*x^3+c*x^2+d*x+e=0 **/.
define @eq4.slv (@eq4.ipt=!charend('/'))
    compute @eq4.b=!@eq4.ipt(1,2)/!@eq4.ipt(1,1).
    compute @eq4.c=!@eq4.ipt(1,3)/!@eq4.ipt(1,1).
    compute @eq4.d=!@eq4.ipt(1,4)/!@eq4.ipt(1,1).
    compute @eq4.e=!@eq4.ipt(1,5)/!@eq4.ipt(1,1).
    @eq3.slv @eq3.ipt={1,-@eq4.c,@eq4.b*@eq4.d-4*@eq4.e,(4*@eq4.c-@eq4.b**2)*@eq4.e-@eq4.d**2}.
    compute @eq4.t1=cmax(@eq.rt).
    compute @eq4.t2=@eq4.b**2+4*@eq4.t1-4*@eq4.c.
    do if (@eq4.t2<=0).
        compute @eq.rtc=0.
        else.
        compute @eq4.tz=0.
        compute @eq4.t3=@eq4.b+sqrt(@eq4.t2).
        compute @eq4.t4=@eq4.t1+(@eq4.b*@eq4.t1-2*@eq4.d)/(sqrt(@eq4.t2)).
        @eq2.slv @eq2.ipt={2,@eq4.t3,@eq4.t4}.
        do if (@eq.rtc=1).
            compute @eq4.tz={@eq4.tz;@eq.rt}.
        end if.
        compute @eq4.t3=@eq4.b-sqrt(@eq4.t2).
        compute @eq4.t4=@eq4.t1-(@eq4.b*@eq4.t1-2*@eq4.d)/(sqrt(@eq4.t2)).
        @eq2.slv @eq2.ipt={2,@eq4.t3,@eq4.t4}.
        do if (@eq.rtc=1).
            compute @eq4.tz={@eq4.tz;@eq.rt}.
        end if.
        compute @eq4.t1=nrow(@eq4.tz).
        do if (@eq4.t1>1).
            compute @eq.rt=@eq4.tz(2:@eq4.t1,1).
            compute @eq.rtc=1.
            else.
            compute @eq.rtc=0.
        end if.
        release @eq4.b,@eq4.c,@eq4.d,@eq4.e,@eq4.t1,@eq4.t2,@eq4.t3,@eq4.t4,@eq4.tz.
    end if.
!enddefine.
/** solve linear, quadratic, cubic, and quartic equation **/.
define @eq.slv (@eq.num=!charend('/'))
    compute @eq.num=!@eq.num.
    compute @eq.pwr=ncol(@eq.num)-1.
    compute @eq.rtc=1.
    do if (@eq.pwr>4 or @eq.pwr<1)./**check errors.
        compute @eq.rtc=0.
        else.
        loop #eq=1 to @eq.pwr.
            do if (@eq.num(1,1)=0).
                do if (@eq.pwr>1).
                    compute @eq.num=@eq.num(1,2:(@eq.pwr+1)).
                    compute @eq.pwr=@eq.pwr-1.
                    else.
                    compute @eq.rtc=0.
                    break.
                end if.
                else.
                break.
            end if.
        end loop.
    end if.
    do if (@eq.rtc=1).
        do if (@eq.pwr=1).
            compute @eq2.rt=-@eq.num(1,2)/@eq.num(1,1).
            else if (@eq.pwr=2).
            @eq2.slv @eq2.ipt=@eq.num.
            else if (@eq.pwr=3).
            @eq3.slv @eq3.ipt=@eq.num.
            else if (@eq.pwr=4).
            @eq4.slv @eq4.ipt=@eq.num.
        end if.
    end if.
    release @eq.num,@eq.pwr.
!enddefine.
/** OLS regression and binary logistic regression **/.
define @reg3 (@rgy=!charend('/') /@rgx=!charend('/') /@rgci=!default(95) !charend('/') 
/@rg1vdr=!default(0) !charend('/') /@rg1hc=!default(9) !charend('/') 
/@rg2cvg=!default(0.001) !charend('/')/@rg2iter=!default(50) !charend('/') 
/@rgALL=!default(1) !charend('/'))
compute rgy=!@rgy. /** dependent variable.
compute rgx=!@rgx. /** independent variables.
compute rgALL=!@rgALL. /** if the value is 0, estimating coefficients only.
compute rgci=!@rgci.
compute rgdichy=ncol(design(rgy)). /** if the value is 2, y is dichotomous.
compute rgnvar=ncol(rgx). /** number of x.
compute rgncase=nrow(rgx). /** number of cases.
compute rgx={make(rgncase,1,1),rgx}. /** add constant to x data.
compute rgym=csum(rgy)/rgncase.  /** mean of y values.
do if (rgdichy=2).
    /** binary logistic regression.
    compute rg2itern=!@rg2iter. /** maximum iterations.
    compute rg2cvg=!@rg2cvg. /** convergence.
    /** recoding for y.
    compute rgymin=cmin(rgy).
    compute rgymax=cmax(rgy).
    compute rgycode=rgy(1,1)-rgymax.
    compute rgy=design(rgy).
    compute rg2y0n=csum(rgy(:,1)).
    do if (rgycode=0).
        compute rgy=rgy(:,1).
        compute rg2y0n=rgncase-rg2y0n.
        else.
        compute rgy=rgy(:,2).
    end if.
    compute rgycode={rgymin,0,rg2y0n;rgymax,1,(rgncase-rg2y0n)}. /** recod info.
    compute rgym=csum(rgy)/rgncase.
    compute rg2p1=make(rgncase,1,0.5).
    compute rgb=make(rgnvar+1,1,0). /** (initial) regression coefficients.
    compute rgb2=rgb. 
    compute rg2LL2=-2*csum(rgy*ln(rgym)+(1-rgy)*ln(1-rgym)). /** (initial) -2LL.
    compute rg2LL1=0.
    compute rg2hstr={0,rg2LL2,t(rgb)}. /** iteration history.
    loop #rg2i=1 to rg2itern.
        compute rgbcov=inv(t(rgx)*mdiag(rg2p1&*(1-rg2p1))*rgx). 
        /** covariance matrix of coefficients.
        compute rgb=rgb2+rgbcov*t(rgx)*(rgy-rg2p1).
        compute rgb2=rgb.
        compute rg2expn=-t(t(rgb2)*t(rgx)).
        compute rg2expmx=cmax(rg2expn).
        do if (rg2expmx>709.77).
            compute rg2exp01=(rg2expn>709.77).
            compute rg2expn=rg2expn&*(1-rg2exp01)+rg2exp01*709.77.
            release rg2exp01.
        end if.
        compute rg2p1=1/(1+exp(rg2expn)).
        release rg2expn.
        compute rg2pmin=cmin(rg2p1).
        compute rg2pmax=cmax(rg2p1).
        do if (rg2pmin=0 or rg2pmax=1).
            loop #cn=1 to rgncase.
                do if (rg2p1(#cn,1)=0).
                    compute rg2p1(#cn,1)=0.0000000001.
                end if.
                do if (rg2p1(#cn,1)=1).
                    compute rg2p1(#cn,1)=0.9999999999.
                end if.
            end loop.
        end if .
        compute rg2LL2=-2*csum(rgy&*ln(rg2p1)+(1-rgy)&*ln(1-rg2p1)).
        compute rg2hstr={rg2hstr;(#rg2i),rg2LL2,t(rgb)}.
        do if (abs(rg2LL1-rg2LL2)<rg2cvg).
            break.
            else.
            compute rg2LL1=rg2LL2.
        end if.
    end loop.
    do if (rgALL=1). /** title.
        compute rgbcov=inv(t(rgx)*mdiag(rg2p1&*(1-rg2p1))*rgx).
        compute rg2se=sqrt(diag(rgbcov)). /** standard errors.
        compute rg2z=rgb/rg2se. /** Z statistics.
        compute rg2pb=2*(1-cdfnorm(abs(rg2z))). /** p values for coefficients.
        compute rg2iterm=nrow(rg2hstr). /** number of iterations.
        compute rg2chisq=rg2hstr(1,2)-rg2LL2. /** Chi-square.
        compute rg2pm=1-chicdf(rg2chisq,rgnvar). /** p value for the model.
        compute rg2CaSR=1-exp(-rg2chisq/rgncase). /** Cox and Snell R-square.
        compute rg2NagR=rg2CaSR/(1-exp(-rg2hstr(1,2)/rgncase)) /** Nagelkerke R-square.
        compute rg2McFR=rg2chisq/rg2hstr(1,2). /** McFadden R-square.
        compute rg2aic=rg2LL2+2*(rgnvar+1).
        compute rg2bic=rg2LL2+rgnvar*ln(rgncase).
        compute rg2aicc=rg2aic+2*(rgnvar+1)*(rgnvar+2)/(rgncase-rgnvar-2).
        compute rg2abic=rg2LL2+rgnvar*ln((rgncase+2)/24).
        compute rgics={rg2aic,rg2bic,rg2aic,rg2abic}.
        /** confidence intervals.
        compute rgci=(100-rgci)/200.
        @ppnd16 idfp=rgci. 
        compute ppnd16=abs(ppnd16).
        compute rgbmx=cmax(rgb2).
        do if (rgbmx>709.77).
            compute rgbexp01=(rgb2>709.77).
            compute rgb2=rgb2&*(1-rgbexp01)+rgbexp01*709.77.
        end if.
        compute rg2coef={rgb,rg2se,rg2z,rg2pb,(rgb-(rg2se*ppnd16)),
            (rgb+(rg2se*ppnd16)),exp(rgb2)}.
        /** standardized coefficients.
        compute rgx={rgy,rg2p1,ln(rg2p1&/(1-rg2p1)),rgx}.
        compute rg2vcv=(sscp(rgx)-sscp(csum(rgx))/rgncase)/(rgncase-1).
        compute rg2sd=sqrt(diag(rg2vcv)).
        compute rg2cr=inv(mdiag(rg2sd(1:2,:)))*rg2vcv(1:2,1:2)*inv(mdiag(rg2sd(1:2,:))).
        compute rg2stdb={rgb&*rg2sd(4:(4+rgnvar),:)*sqrt(3)/(4*artan(1)),
            rgb&*rg2sd(4:(4+rgnvar),:)*rg2cr(1,2)/rg2sd(3,1)}. 
        compute rg2coef={rg2coef,rg2stdb}.
        compute rg2ifm={rg2cvg,(rg2iterm-1),rg2hstr(rg2iterm,2),
            rg2chisq,rgnvar,rg2pm,rg2CaSR,rg2NagR,rg2McFR}.
        /** Classification.
        compute rg2predg=rnd(rg2p1). 
        compute rg2g4=rgy+rg2predg/2.
        compute rg2g4={rg2g4;0;0.5;1;1.5}.
        @sortv sortv=rg2g4/sorttp=1.
        compute rg2g4=csum(design(rg2g4)).
        compute rg2class={rg2g4(1,1:2);rg2g4(1,3:4)}. 
        compute rg2class=rg2class-1.
        compute rg2class={rg2class,rsum(rg2class)}.
        compute rg2class={rg2class;csum(rg2class)}.
        do if (rg2class(1,3)<>0).
            compute rg2class(1,3)=rg2class(1,1)/rg2class(1,3)*100.
        end if.
        do if (rg2class(2,3)<>0).
            compute rg2class(2,3)=rg2class(2,2)/rg2class(2,3)*100.
        end if.
        compute rg2class(3,3)=(rg2class(1,1)+rg2class(2,2))/rg2class(3,3)*100.
        do if (rg2class(3,1)<>0).
            compute rg2class(3,1)=rg2class(1,1)/rg2class(3,1)*100.
        end if.
        do if (rg2class(3,2)<>0).
            compute rg2class(3,2)=rg2class(2,2)/rg2class(3,2)*100.
        end if.
        compute rgycodep={'Raw','Analysis','N'}.
        compute rgicsp={'AIC','BIC','AICc','aBIC'}.
        compute rgbp={'Coeff','S.E.','Z','p','LLCI','ULCI','OR','Std.C.1','Std.C.2'}.
        compute rg2ifmp={'Convergence','N_iterat','-2LL','Chi-sq','df',
            'p','CoxSnell','Nagelkerke','McFadden'}.
        compute rg2clspc={'Pred.0','Pred.1','%'}.
        compute rg2clspr={'Raw.0','Raw.1','%','       .'}.
    end if.
else.
    compute rg1vdr=!@rg1vdr.
    compute rg1hc=!@rg1hc.
    compute rg1ixtx=inv(t(rgx)*rgx).
    compute rgb=rg1ixtx*t(rgx)*rgy. /** regression coefficients.
    do if (rgALL=1).
        compute rg1predy=rgx*rgb. /** prdicted y values.
        compute rg1yss=cssq(rgy-rgym).  /** y - sum of squares.
        compute rg1e=rgy-rg1predy. /** residuals.
        compute rg1ess=cssq(rg1e). /** error - sum of squares.
        compute rg1rss=rg1yss-rg1ess. /** regression - sum of squares.
        compute rg1df1=rgnvar.  
        compute rg1df2=rgncase-rgnvar-1.
        compute rg1rssm=rg1rss/rg1df1. 
        compute rg1essm=rg1ess/rg1df2.
        compute rg1f=rg1rssm/rg1essm. /** f-ratio.
        compute rg1pm=1-fcdf(rg1f,rg1df1,rg1df2). /** p value for the model.
        compute rg1rsq=rg1rss/rg1yss. /** R-square.
        compute rg1arsq=1-(1-rg1rsq)*(rgncase-1)/rg1df2. /** adjusted R-square.
        compute rg1ifm={rg1rsq,rg1arsq,rg1essm,rg1f,rg1df1,rg1df2,rg1pm}.
        compute rg1aic=rgncase*ln(rg1ess/rgncase)+2*(rgnvar+1). 
        compute rg1bic=rgncase*ln(rg1ess/rgncase)+(rgnvar+1)*ln(rgncase).
        compute rg1aicc=rg1aic+2*(rgnvar+1)*(rgnvar+2)/(rgncase-rgnvar-2).
        compute rg1abic=rgncase*ln(rg1ess/rgncase)+ln((rgncase+2)/24)*(rgnvar+1). 
        compute rgics={rg1aic,rg1bic,rg1aicc,rg1abic}. 
        compute rgbcov=rg1ixtx*rg1essm. 
        /** covariance matrix of coefficients.
        compute rg1se0=sqrt(diag(rgbcov)).
        do if (!@rg1hc<9). /** robust standard errors.
            compute rg1e2=rg1e&**2.
            do if (!@rg1hc=2).
                compute rg1e2=rg1e2&/(make(ncase,1,1)-diag(rgx*rg1ixtx*t(rgx))).
                else if (!@rg1hc=3).
                compute rg1e2=rg1e2&/((make(ncase,1,1)-diag(rgx*rg1ixtx*t(rgx)))&**2).
                else if (!@rg1hc=4).
                compute rg1hat=rgx*rg1ixtx*t(rgx).
                compute rg1tmp1={make(rgncase,1,4),(rgncase*diag(rg1hat)/(rgnvar+1))}.
                compute rg1e2=rg1e2&/((make(ncase,1,1)-diag(rg1hat))&**rmin(rg1tmp1)).
            end if.
            compute rgbcov=rg1ixtx*t(rgx)*mdiag(rg1e2)*rgx*rg1ixtx.
            do if (!@rg1hc=1).
                compute rgbcov=rgbcov*rgncase/(rg1df2).
            end if.
        end if.
        compute rg1se=sqrt(diag(rgbcov)). /** standard errors.
        compute rg1t=rgb/rg1se. /** t values for coefficients.
        compute rg1pb=2*(1-tcdf(abs(rg1t),rg1df2)). /** p values for coefficients.
        /** standardized coefficients.
        compute rgx(:,1)=rgy.
        compute rg1vcv=1/(rgncase-1)*(sscp(rgx)-((t(csum(rgx))*csum(rgx))/rgncase))).
        compute rg1sd=sqrt(diag(rg1vcv)).
        compute rg1stdb=rgb&*rg1sd/(rg1sd(1,1)).
        compute rg1stdb(1,1)=0.
        @qtidf qtdf=rg1df2/qtprcn=rgci.
        compute rg1coef={rgb,rg1se,rg1t,rg1pb,(rgb-rg1se*qtout),(rgb+rg1se*qtout),rg1stdb}.
        do if (!@rg1vdr=1).
            do if (rgnvar=1). /** VIF.
                compute rg1vif={0;0}.
                else.
                compute rg1x=rgx(:,2:(rgnvar+1)).
                compute rg1vcv=1/(rgncase-1)*(sscp(rg1x)-((t(csum(rg1x))*csum(rg1x))/rgncase))).
                compute rg1d=inv(mdiag(sqrt(diag(rg1vcv)))).
                compute rg1cr=rg1d*rg1vcv*rg1d.
                compute rg1vif={0;diag(inv(rg1cr)*rg1cr*inv(rg1cr))}.
                release rg1x.
            end if.
            compute rg1rsqc=((rgb&/rg1se0)*((1-rg1rsq)/rg1df2)**0.5)&**2.
            compute rg1rsqc={0;rg1rsqc(2:(rgnvar+1),1)}. /** Delta R-square.
            compute rg1coef={rg1coef,rg1vif,rg1rsqc}.
        end if.
        compute rgicsp={'AIC','BIC','AICc','aBIC'}.
        compute rgbp={'Coeff','S.E.','t','p','LLCI','ULCI','Std.Coeff','VIF','DeltaRsq'}.
        compute rg1ifmp={'R-sq','Adj.R-sq','MSE','F','df1','df2','p'}.
    end if.
end if.
release rgy,rgx.
!enddefine. 
/** inverse cdf of t distribution **/.
define @qtidf(qtdf=!charend('/') / qtprcn=!charend('/'))
compute qtn=!qtdf.
compute qtp=1-(!qtprcn/100).
do if (qtn<1).
    compute qtn=1.
end if .
do if (qtp>=1 or qtp<=0).
    compute qtp=0.05.
end if.
do if (qtn=2).
    compute qtout=sqrt(2/(qtp*(2-qtp))-2).
    ELSE.
    compute qthpi=2*artan(1).
    do if (qtn=1).
        compute qtp=qtp*qthpi.
        compute qtout=cos(qtp)/sin(qtp).
        ELSE.
        compute qta=1/(qtn-0.5).
        compute qtb=48/(qta**2).
        compute qtc=((20700*qta/qtb-98)*qta-16)*qta+96.36.
        compute qtd=((94.5/(qtb+qtc)-3)/qtb+1)*sqrt(qta*qthpi)*qtn.
        compute qtx=qtd*qtp.
        compute qty=qtx**(2/qtn).
        do if(qty>(0.05+qta)).
            compute qtx=qtp*0.5.
            @ppnd16 idfp=qtx.
            compute qtx=ppnd16.
            compute qty=qtx**2.
            do if (qtn<5).
                compute qtc=qtc+0.3*(qtn-4.5)*(qtx+0.6).
            end if.
            compute qtc=(((0.05*qtd*qtx-5)*qtx-7)*qtx-2)*qtx+qtb+qtc.
            compute qty=(((((0.4*qty+6.3)*qty+36)*qty+94.5)/qtc-qty-3)/qtb+1)*qtx.
            compute qty=qta*(qty**2).
            do if (qty>0.01).
                do if (qty>709.77).
                    compute qty=709.77.
                end if.
                compute qty=exp(qty)-1.
                ELSE. 
                compute qty=((qty+4)*qty+12)*qty*qty/24+qty.
            end if.
        ELSE.
            compute qty=((1/(((qtn+6)/(qtn*qty)-0.089*qtd-0.822)*(qtn+2)*3)+0.5
                /(qtn+4))*qty-1)*(qtn+1)/(qtn+2)+1/qty.
        end if.
        compute qtout=sqrt(qtn*qty).
    end if.
end if.
!enddefine.
/** PCR **/.
define @PCAreg2(pcr_y=!charend('/') /pcr_x = !charend('/')/pcr_n = !charend('/'))
compute dataall={!pcr_y,!pcr_x}.
compute ncase=nrow(dataall).
compute nvars=ncol(dataall).
compute nfac=!pcr_n.
compute nvarx=nvars-1.
compute xysdm=make(2,nvars,-999).
compute stdall=make(ncase,nvars,-999).
compute tdata1=1/(ncase-1).
compute tdata2=tdata1*(sscp(dataall)-((t(csum(dataall))*csum(dataall))/ncase)).
compute tdata3=inv(mdiag(sqrt(diag(tdata2)))).
compute corrall=tdata3*tdata2*tdata3.
compute tdata1=sqrt(t(diag(tdata2))).
compute tdata2=csum(dataall)/ncase.
compute xysdm={tdata2;tdata1}.
loop #i=1 to nvars.
   compute stdall(:,#i)=(dataall(:,#i)-tdata2(1,#i))/tdata1(1,#i).
end loop.
compute stdx=stdall(:,2:nvars).
compute stdy=stdall(:,1).
compute rawx=dataall(:,2:nvars).
compute rawy=dataall(:,1).
compute corrx=corrall(2:nvars,2:nvars).
call eigen(corrx,evector,evalue).
compute cpmatrix=evector*sqrt(mdiag(evalue)).
compute scmatrix=make(nvarx,nvarx,0).
loop #i=1 to nvarx.
    compute scmatrix(:,#i)=cpmatrix(:,#i)/evalue(#i,1).
end loop.
compute cpscore=make(ncase,nvarx,0).
/** regression factor score **/. 
loop #i=1 to nvarx.
    loop #i2=1 to nvarx.
        compute cpscore(:,#i)=cpscore(:,#i)+stdx(:,#i2)*scmatrix(#i2,#i). 
        do if (#i2=nvarx).
            compute cpscore(:,#i)=cpscore(:,#i)*sqrt(evalue(#i,1)).
        end if.
    end loop.
end loop.
@Polsb ols_y=stdy/ols_x=cpscore(:,1:nfac).
compute stdxcoef=make(nvarx,1,0).
compute xcoef=stdxcoef.
compute rawconst=xysdm(1,1).
loop #i=1 to nvarx.
    loop #i2=1 to nfac.
        compute stdxcoef(#i,1)=stdxcoef(#i,1)+evector(#i,#i2)*olsb(#i2+1,1).
        do if (#i2=nfac).
            compute xcoef(#i,1)=stdxcoef(#i,1)/xysdm(2,#i+1)*xysdm(2,1).
        end if.
    end loop.
    do if (#i=nvarx).
        loop #i3=1 to nvarx.
            compute rawconst=rawconst-xcoef(#i3,1)*xysdm(1,#i3+1).
        end loop.
    end if.
end loop.
compute ssy=csum((rawy-xysdm(1,1))&**2).
loop #i=1 to nvarx.
    compute rawy=rawy-rawx(:,#i)*xcoef(#i,1).
    do if (#i=nvarx).
        compute rawy=rawy-rawconst.
    end if.
end loop.
compute sse=csum((rawy-csum(rawy)/ncase)&**2).
compute ssm=ssy-sse.
compute modelr2=ssm/ssy.
compute ftest=make(3,5,0).
compute ftest(1:3,1)={ssm;sse;ssy}.
compute ftest(1:3,2)={nvarx;(ncase-1-nvarx);(ncase-1)}.
compute ftest(1,3)=ftest(1,1)/ftest(1,2).
compute ftest(2,3)=ftest(2,1)/ftest(2,2).
compute ftest(1,4)=ftest(1,3)/ftest(2,3).
compute ftest(1,5)=1-fcdf(ftest(1,4),nvarx,ftest(2,2)).
compute xcoef={rawconst;xcoef}.
!enddefine.
/** inverse cdf of normal distribution **/.
define @ppnd16(idfp=!charend('/'))
compute invnp=!idfp.
compute invna0=3.3871328727963666080.
compute invna1=133.14166789178437745.
compute invna2=1971.5909503065514427.
compute invna3=13731.693765509461125.
compute invna4=45921.953931549871457.
compute invna5=67265.770927008700853.
compute invna6=33430.575583588128105.
compute invna7=2509.0809287301226727.
compute invnb1=42.313330701600911252.
compute invnb2=687.18700749205790830.
compute invnb3=5394.1960214247511077.
compute invnb4=21213.794301586595867.
compute invnb5=39307.895800092710610.
compute invnb6=28729.085735721942674.
compute invnb7=5226.4952788528545610.
compute invnc0=1.42343711074968357734.
compute invnc1=4.63033784615654529590.
compute invnc2=5.76949722146069140550.
compute invnc3=3.64784832476326738258.
compute invnc4=1.27045825245236838258.
compute invnc5=0.241780725177450611770.
compute invnc6=0.0227238449892691845833.
compute invnc7=0.000774545014278341407640.
compute invnd1=2.05319162663775882187.
compute invnd2=1.67638483018380384940.
compute invnd3=0.689767334985100004550.
compute invnd4=0.148103976427480074590.
compute invnd5=0.0151986665636164571966.
compute invnd6=0.000547593808499534494600.
compute invnd7=0.105075007164441684324.
compute invne0=6.65790464350110377720.
compute invne1=5.46378491116411436990.
compute invne2=1.78482653991729133580.
compute invne3=0.296560571828504891230.
compute invne4=0.0265321895265761230930.
compute invne5=0.00124266094738807843860.
compute invne6=0.0000271155556874348757815.
compute invne7=0.000000201033439929228813265.
compute invnf1=0.599832206555887937690.
compute invnf2=0.136929880922735805310.
compute invnf3=0.0148753612908506148525.
compute invnf4=0.000786869131145613259100.
compute invnf5=0.0000184631831751005468180.
compute invnf6=0.000000142151175831644588870.
compute invnf7=0.204426310338993978564.
compute invnq=invnp-0.5.
do if (abs(invnq) <= 0.425).
    compute invnr=0.180625-invnq**2.
    compute ppnd16fz=invnq*(((((((invna7*invnr+invna6)*invnr+invna5)*invnr+
        invna4)*invnr+invna3)*invnr+invna2)*invnr+invna1)*invnr+invna0).
    compute ppnd16fm=(((((((invnb7*invnr+invnb6)*invnr+invnb5)*invnr+invnb4)*
        invnr+invnb3)*invnr+invnb2)*invnr+invnb1)*invnr+1).
    compute ppnd16=ppnd16fz/ppnd16fm.
    else.
    do if (invnq<=0).
        compute invnr=invnp.
        ELSE.  
        compute invnr=1-invnp.
    end if.
    do if (invnr<=0).
        compute ppnd16=0.
        else.
        compute invnr=sqrt(-ln(invnr)).
        do if (invnr<=5).
            compute invnr=invnr-1.6.
            compute ppnd16fz=(((((((invnc7*invnr+invnc6)*invnr+invnc5)*invnr+
                invnc4)*invnr+invnc3)*invnr+invnc2)*invnr+invnc1)*invnr+invnc0).
            compute ppnd16fm=(((((((invnd7*invnr*0.00000001+invnd6)*invnr+
                invnd5)*invnr+invnd4)*invnr+invnd3)*invnr+invnd2)*invnr+invnd1)*invnr+1).
            compute ppnd16=ppnd16fz/ppnd16fm.
            else.
            compute invnr=invnr-5. 
            compute ppnd16fz=(((((((invne7*invnr+invne6)*invnr+invne5)*invnr+
                invne4)*invnr+invne3)*invnr+invne2)*invnr+invne1)*invnr+invne0).
            compute ppnd16fm=(((((((invnf7*invnr*0.00000000000001+invnf6)*
                invnr+invnf5)*invnr+invnf4)*invnr+invnf3)*invnr+invnf2)*invnr+invnf1)*invnr+1).
            compute ppnd16=ppnd16fz/ppnd16fm.
        end if.
        do if (invnq<0).
            compute ppnd16=-ppnd16.
        end if.
    end if.
end if.
!enddefine.
/** Bootstrap's output **/.
define @btsot(btrawd=!charend('/')/btalld=!charend('/')/btcitp=!default(1) !charend('/')
    /btnewnm=!charend('/')/btcih=!default(97.5) !charend('/')/btcil=!default(2.5) !charend('/'))
    compute btcinvar=ncol(!btrawd).
    compute btcin=nrow(!btalld).
    compute !btnewnm=make(6,btcinvar,-999).
    compute !btnewnm(1,:)=!btrawd.
    compute !btnewnm(2,:)=csum(!btalld)/btcin.
    loop #btc=1 to btcinvar.
        compute btmpdt=!btalld(:,#btc).
        compute btmpdt2=make(btcin,1,!btnewnm(2,#btc)).
        compute !btnewnm(4,#btc)=(csum((btmpdt-btmpdt2)&**2)/(btcin-1))**0.5.
        @prcn prcny=btmpdt/prcnp=50.
        compute !btnewnm(3,#btc)=percnv.
        do if (!btcitp=1).
            @prcn prcny=btmpdt/prcnp=!btcil.
            compute !btnewnm(5,#btc)=percnv. 
            @prcn prcny=btmpdt/prcnp=!btcih.
            compute !btnewnm(6,#btc)=percnv.
            else.    
            @bcci bcciy=btmpdt/bccih=bcrcih/bcciraw=!btrawd(1,#btc).
            compute !btnewnm(5,#btc)=bccilp. 
            compute !btnewnm(6,#btc)=bccilh. 
        end if.
    end loop.
!enddefine.
/** Bias-corrected CI **/.
define @bcci(bcciy=!charend('/')/bccih=!default(97.5) !charend('/')/bcciraw=!charend('/'))
compute bccixx=!bcciy.
compute bccip3=!bccih.
@ppnd16 idfp=bccip3.
compute bccirz=ppnd16.
compute bccin=nrow(bccixx).
compute bcciul=99.94-98/bccin.
compute bccill=0.01+99/bccin.
compute bccinc=0.
loop #bc=1 to bccin.
    do if (bccixx(#bc,1)<!bcciraw).
        compute bccinc=bccinc+1.
    end if.
end loop.
compute bccinp=bccinc/(bccin+1).
@ppnd16 idfp=bccinp.
compute bccimz=ppnd16.
compute bccilz=2*bccimz-bccirz.
compute bccihz=2*bccimz+bccirz.
compute bccilp=cdfnorm(bccilz)*100.
compute bccihp=cdfnorm(bccihz)*100.
do if (bccihp>bcciul).
    compute bccihp=bcciul.
    do if (csum(ero)<>12).
        compute ero={ero;12}.
    end if.
end if.
do if (bccihp<bccill).
    compute bccihp=bccill.
    do if (csum(ero)<>12).
        compute ero={ero;12}.
    end if.
end if.
do if (bccilp<bccill).
    compute bccilp=bccill.
    do if (csum(ero)<>12).
        compute ero={ero;12}.
    end if.
end if.
do if (bccilp>bcciul).
    compute bccilp=bcciul.
    do if (csum(ero)<>12).
        compute ero={ero;12}.
    end if.
end if.
@prcn prcny=bccixx/prcnp=bccilp.
compute bccilp=percnv.
@prcn prcny=bccixx/prcnp=bccihp.
compute bccilh=percnv.
!enddefine.
/** Data sorting **/.
define @sortv(sortv=!charend('/')/sorttp=!default(1) !charend('/'))
    compute sortdt=!sortv.
    compute sorttp=!sorttp.    
    compute sort2=grade(sorttp&*sortdt).
    compute sortdt(sort2)=sortdt.
    compute !sortv=sortdt.
    release sortdt.
!enddefine.
/** Percentile **/.
define @prcn(prcny=!charend('/')/prcnp=!charend('/'))
compute percnp=!prcnp.
compute percny=!prcny.
compute percnn=nrow(percny).
compute percnnn=(percnn+1)*percnp/100.
compute percni=trunc(percnnn).
compute percni2=percni+1.
compute percnd=percnnn-percni.
do if (percnd =0).
    compute  percnv=percny(percni,1).
    else.
    compute  percnv=percnd*percny(percni2,1)+(1-percnd)*percny(percni,1).
end if.
!enddefine.
/** Check constant and correlation **/.
define @candr(crdata=!charend('/'))
    compute crdata=!crdata.
    compute crnvar=ncol(crdata).
    compute candrn=0.
    loop #i=1 to crnvar.
        compute tdata=crdata(1:ncase,#i).
        compute tmax=cmax(tdata).
        compute tmin=cmin(tdata).
        do if (tmax=tmin).
            compute candrn=1.
        end if.
    end loop.
    do if (candrn=0).
        compute tpar1=1/(ncase-1).
        compute corvcv=tpar1*(sscp(crdata)-((t(csum(crdata))*csum(crdata))/ncase)).
        compute cord=inv(mdiag(sqrt(diag(corvcv)))).
        compute corcr=cord*corvcv*cord.
        loop #ckcr1=1 to crnvar.
            loop #ckcr2 =1 to crnvar.
                do if (#ckcr1 <> #ckcr2).
                    do if (corcr(#ckcr1,#ckcr2)=1).
                        compute candrn=2.
                    end if.
                end if.
            end loop.
        end loop.
    end if.
    do if (candrn=0).
        compute ckeval=eval(corcr).
        do if (cmin(ckeval)<0.00001).
            compute candrn=3.
        end if.
    end if.
    release crdata,crnvar,tdata,tmax,tmin,tpar1,corvcv,cord,corcr,ckeval.
!enddefine.
/** OLS-B 4 PCA **/.
define @Polsb(ols_y=!charend('/') /ols_x = !charend('/')).
    compute olsy=!ols_y.
    compute ncase=nrow(olsy).
    compute olsx={make(ncase,1,1),!ols_x}.
    compute olsb=ginv(olsx)*olsy.
!enddefine.
define ereg1.0.1 (emod=!default(1) !charend('/')/e_y=!charend('/')/e_x=!charend('/')
/seed=!default(random) !charend('/')/ci=!default(95) !charend('/')
/pcrtype=!default(1) !charend('/')/pcrnum=!default(1) !charend('/')
/pcrdes=!default(0) !charend('/')/pcrcor=!default(0) !charend('/')
/pcrallrb=!default(0) !charend('/')/pcreva=!default(0) !charend('/')
/pcreve=!default(0) !charend('/')/pcrcm=!default(0) !charend('/')
/pcrcsm=!default(0) !charend('/')/pcrplot1=!default(0) !charend('/')
/pcrplot2=!default(0) !charend('/')/pcrbtn=!default(0) !charend('/')
/pcrbtm=!default(0) !charend('/')/dcms=!default(f12.8) !charend('/')
/pcrplot3=!default(0) !charend('/')/rrkmin=!default(0) !charend('/')
/rrkmax=!default(1) !charend('/')/rrkinc=!default(0) !charend('/')
/rrkstd=!default(0) !charend('/')/rrkvif=!default(0) !charend('/')
/rrkfit=!default(0) !charend('/')/rrkm=!default(1) !charend('/')
/rrkfin=!default(0) !charend('/')/rrbtn=!default(1000) !charend('/')
/rrbtm=!default(1) !charend('/')/rrplot1=!default(0) !charend('/')
/rrplot2=!default(0) !charend('/')/rrplot3=!default(0) !charend('/')
/rrdes=!default(0) !charend('/')/rrcor=!default(0) !charend('/')
/rgpm=!default(1) !charend('/')/rgpdn=!default(11) !charend('/')
/mrres=!default(0) !charend('/')/mrdsts=!default(0) !charend('/')
/mrdcasen=!default(5) !charend('/')/mrcoefev=!default(0) !charend('/')
/mrca=!default(0) !charend('/')/mrda=!default(0) !charend('/')
/mrallsub=!default(0) !charend('/')/mrassot=!default(1) !charend('/')
/mrdes=!default(0) !charend('/')/mrcovb=!default(0) !charend('/')
/mrheter=!default(0) !charend('/')/mrrbste=!default(9) !charend('/')
/mrplot1=!default(0) !charend('/')/mrplot2=!default(0) !charend('/')
/qedttp=!default(1) !charend('/')/qejn=!default(0) !charend('/')
/qerbst=!default(9) !charend('/')/qepap=!default(1) !charend('/')
/qejnp=!default(0) !charend('/')/qepapv=!default(0) !charend('/')
/spltp=!default(1) !charend('/')/splktp=!default(1) !charend('/')
/splqv=!default(0) !charend('/')/splrbst=!default(9) !charend('/')
/splsave=!default(0) !charend('/')/qesave=!default(0) !charend('/')
/rgpmol=!default(0) !charend('/')/rgpmov=!default(0) !charend('/')
/rgpemp=!default(1) !charend('/')/rgpemq=!default(0) !charend('/')
/rgprbst=!default(9) !charend('/')/splplot=!default(0) !charend('/')
/splplotn=!default(0) !charend('/')/splemmp=!default(0) !charend('/')
/modm=!default(6) !charend('/')/moddttp=!default(1) !charend('/')
/modsave=!default(1) !charend('/')/modrbst=!default(9) !charend('/')
/modpap=!default(1) !charend('/')/modpapv=!default(0) !charend('/')
/modjn=!default(0) !charend('/')/modjnp=!default(0) !charend('/')
/modjn2=!default(0) !charend('/')/modjnx=!default(0) !charend('/')
/lgconver=!default(0.001) !charend('/')/modjnn=!default(20) !charend('/'))
set seed=!seed.
compute EReg_ID=$casenum.
execute.
MATRIX.
!let !tt=('*********************************************'+
    '********************************************')
!let !tt=!quote(!tt)
print /title=!tt/format=a8.
print /title="EReg version 1.0.1"/foramt=a8/space=0.
print /title='Copyright (C) 2020-2021 by Qiu Zongman'/format=a8.
print /title=!tt/format=a8.
/** extract data and settings **/.
get raw1/file=*/var=!e_y !e_x /name=nmall/missing=omit.
get raw2/file=*/var=!e_y /name=nmall2/missing=omit.
get eregid/file=*/vars=EReg_ID !e_y !e_x /missing=omit.
compute eregid=eregid(:,1).
compute nvarall=ncol(raw1).
compute nvaryw=ncol(raw2).
compute ncase=nrow(raw1).
compute lgconver=!lgconver.
do if (ncase<>1).
    compute ncase2=1/(ncase-1).
end if.
compute emod=!emod.
compute edty=raw1(:,1).
compute dichy=ncol(design(edty)).
compute edtx1=raw1(:,2).
compute edtx1max=cmax(edtx1).
compute edtx1min=cmin(edtx1).
compute edtx2=raw1(:,2:nvarall).
compute nvarx2=nvarall-1.
compute nmy=nmall(1,1).
compute nmx1=nmall(1,2).
compute nmx2=nmall(1,2:nvarall).
compute nmx2={nmx2,'       .'}.
compute nvarco=nvarall-2.
do if (nvarco>0).
    compute edtco=raw1(:,3:nvarall).
    compute nmco=nmall(1,3:nvarall).
end if.
do if (!quote(!seed)="random").
   compute seed={"random"}.
   else.
   compute seed={!quote(!seed)}.
end if.
compute eregci=!ci.
do if (eregci<50).
   compute eregci=95.
end if.
do if (eregci>=100).
   compute eregci=95.
end if.
compute eregcih=50+eregci/2.
compute eregcil=100-eregcih.
compute bcrcih=eregcih/100.
compute ero=0.
@candr crdata=raw1.
compute ero={ero;candrn}.
compute sero=csum(ero).
compute btero=0.
compute pcrbtm=!pcrbtm.
compute rrkcrt=0.
compute prtrslt={'Results:'}.
compute colnm={"COL1","COL2","COL3","COL4","COL5","COL6","COL7","COL8","COL9","COL10",
    "COL11","COL12","COL13","COL14","COL15","COL16","COL17","COL18","COL19","COL20","COL21",
    "COL22","COL23","COL24","COL25","COL26","COL27","COL28","COL29","COL30","COL31","COL32",
    "COL33","COL34","COL35","COL36","COL37","COL38","COL39","COL40","COL41","COL42","COL43",
    "COL44","COL45","COL46","COL47","COL48","COL49","COL50","COL51","COL52","COL53","COL54",
    "COL55","COL56","COL57","COL58","COL59","COL60","COL61","COL62","COL63","COL64","COL65",
    "COL66","COL67","COL68","COL69","COL70","COL71","COL72","COL73","COL74","COL75","COL76",
    "COL77","COL78","COL79","COL80","COL81","COL82","COL83","COL84","COL85","COL86","COL87",
    "COL88","COL89","COL90","COL91","COL92","COL93","COL94","COL95","COL96","COL97","COL98",
    "COL99","COL100","COL101","COL102","COL103","COL104","COL105","COL106","COL107","COL108",
    "COL109","COL110","COL111","COL112","COL113","COL114","COL115","COL116","COL117","COL118",
    "COL119","COL120","COL121","COL122","COL123","COL124","COL125","COL126","COL127","COL128",
    "COL129","COL130","COL131","COL132","COL133","COL134","COL135","COL136","COL137","COL138",
    "COL139","COL140","COL141","COL142","COL143","COL144","COL145","COL146","COL147","COL148",
    "COL149","COL150"}.
do if (emod=6).
        print /title='Module :  Principal Component Regression (PCR) '/format=a8/space=0.
    else if (emod=2).
        print /title='Module :  Multiple Regression'/format=a8/space=0.
    else if (emod=5).
        print /title='Module :  Ridge Regression (RR) '/format=a8/space=0.
    else if (emod=8).
        print /title='Module :  Spline Regression'/format=a8/space=0.
    else if (emod=3).
        print /title='Module :  Quadratic Effect'/format=a8/space=0.
    else if (emod=1).
        print /title='Module :  Visualization'/format=a8/space=0.
    else if (emod=4).
        print /title='Module :  Moderation'/format=a8/space=0.
end if.
!let !tt = !concat('Dependent Variable (DV): ',!e_y).
!let !tt = !quote(!tt).
print /title=!tt/format=a8/space=1.
!let !tt = !concat('Independent Variables (IVs): ',!e_x).
!let !tt = !quote(!tt).
print ncase/title=!tt/rname={'N      :'}/foramt=f8.0/space=0.
!let !tt=!concat('Confidence Intervals (CIs) :',!ci,'%')
!let !tt=!quote(!tt)
print /title=!tt/space=0.
!let !rgline=!concat('************************************ Regression '+
    'Model ***********************************')
!let !rgline=!quote(!rgline)
!let !rgytitle=!concat('====================> Outcome Variable : ',!e_y)
!let !rgytitle=!quote(!rgytitle)
!let !rgycode=!concat('>>>>> Recoding of ',!e_y)
!let !rgycode=!quote(!rgycode)
!let !plotitle=!concat('************************************* Code '+
    'for Plot *************************************')
!let !plotitle=!quote(!plotitle)
!let !plonote=!concat('Note   : Paste the text below into an SPSS syntax ','editor and run.')
!let !plonote=!quote(!plonote)
do if (emod=1). /** Visualization.
    do if (dichy<>2).
        compute rgprbst=!rgprbst.
        do if (rgprbst<>9).
            !let !tt = !concat('Standard errors for regression coefficients : HC',!rgprbst)
            !let !tt = !quote (!tt)
            print /title=!tt/format=a8.
        end if .
    end if.
    print /title 'Covariates appearing in the model are evaluated at mean value'/format=a8/space=0.
    compute tdata1=nmx2(1,1).
    print tdata1/title 'Focal IV (X) is '/format=a8.
    compute rgpm=!rgpm.
    do if (rgpm>5.5).
        do if (nvarx2=1).
            print /title 'Moderator does not exist. '/format=a8/space=0.
            compute ero={ero;14}.
        else.
            compute tdata2={nmx2(1,2)}.
            print tdata2/title 'Moderator (M) is '/format=a8/space=0.
            compute intnm={'X^2','=',tdata1,'*',tdata1,' ',' ';'M^2','=',tdata2,'*',tdata2,' ',' ';
                'X*M','=',tdata1,'*',tdata2,' ',' ';'X^2*M','=',tdata1,'*',
                tdata1,'*',tdata2;'X*M^2','=',tdata1,'*',tdata2,'*',tdata2}.
            do if (rgpm=6).
                compute intnmp=intnm(3,:).
                else if (rgpm=7).
                compute intnmp={intnm(3,:);intnm(1,:)}.
                else if (rgpm=8).
                compute intnmp={intnm(3,:);intnm(1,:);intnm(4,:)}.
                else if (rgpm=9).
                compute intnmp={intnm(3,:);intnm(2,:)}.
                else if (rgpm>=10).
                compute intnmp={intnm(3,:);intnm(2,:);intnm(5,:)}.
            end if.
            print intnmp/title 'Product term(s) key:'/rname={'List   :'}/format=a8/space=0.
        end if.
    end if.
    compute rgpdn=!rgpdn.
    compute rgpxmax=cmax(edtx2).
    compute rgpxmin=cmin(edtx2).
    compute sero=csum(ero).
    do if (sero=0).
        do if (rgpm<5.5).
            compute tdata1={nmx1,'       :'}.
            else.
            compute tdata1={nmx2(1,1:2),'       :'}.
            compute tdata2={t(rgpxmin(1,1:2)),t(rgpxmax(1,1:2))}.
        end if.
    end if.
    do if (rgpm<5.5). /** check errors - sample size.
        compute tdata1=nvarx2+rgpm-1.
    else if (rgpm=6).
        compute tdata1=nvarx2+1.
    else if (rgpm=7).
        compute tdata1=nvarx2+2.
    else if (rgpm=8).
        compute tdata1=nvarx2+3.
    else if (rgpm=9).
        compute tdata1=nvarx2+2.
    else if (rgpm=10).
        compute tdata1=nvarx2+3.
    end if.
    do if (tdata1>(ncase-2)).
        compute ero={ero;7}.
    end if.
    compute sero=csum(ero).
    do if (rgpdn<11).
       compute rgpdn=11.
    end if.
    do if (sero=0).
        do if (rgpm<5.5).
            compute tdata1={edtx1,edtx1&**2,edtx1&**3,edtx1&**4,edtx1&**5}.
            compute tdata1=tdata1(:,1:rgpm).
            compute tdatanm={'Constant',nmx1,'X^2','X^3','X^4','X^5'}.
            compute tdatanm=tdatanm(1,1:(rgpm+1)).
            do if (nvarx2>1).
                compute tdata1={tdata1,edtco}.
                compute tdatanm={tdatanm,nmx2(1,2:nvarx2)}.
            end if.
        else.
            compute tdata1=edtx2(:,2).
            do if (rgpm=6).
                compute tdata1={edtx1,tdata1,edtx1&*tdata1}.
                compute tdatanm={'Constant',nmx1,nmx2(1,2),'X*M'}.
            else if (rgpm=7).
                compute tdata1={edtx1,edtx1&**2,tdata1,edtx1&*tdata1}.
                compute tdatanm={'Constant',nmx1,'X^2',nmx2(1,2),'X*M'}.
            else if (rgpm=8).
                compute tdata1={edtx1,edtx1&**2,tdata1,edtx1&*tdata1,edtx1&**2&*tdata1}.
                compute tdatanm={'Constant',nmx1,'X^2',nmx2(1,2),'X*M','X^2*M'}.
            else if (rgpm=9).
                compute tdata1={edtx1,tdata1,tdata1&**2,edtx1&*tdata1}.
                compute tdatanm={'Constant',nmx1,nmx2(1,2),'M^2','X*M'}.
            else if (rgpm=10).
                compute tdata1={edtx1,tdata1,tdata1&**2,edtx1&*tdata1,edtx1&*tdata1&**2}.
                compute tdatanm={'Constant',nmx1,nmx2(1,2),'M^2','X*M','X*M^2'}.
            end if.
            do if (nvarx2>2).
                compute tdata1={tdata1,edtx2(:,3:nvarx2)}.
                compute tdatanm={tdatanm,nmx2(1,3:nvarx2)}.
            end if.
        end if. 
        @candr crdata={tdata1,edty}.
        do if (candrn>0).
            compute ero={ero;3}.
            else.
            compute tdatax={edty,edtx2}.
            compute rgpvcv=ncase2*(sscp(tdatax)-((t(csum(tdatax))*csum(tdatax))/ncase)).
            compute rgpsd=sqrt(diag(rgpvcv)).
            compute rgpmean=csum(tdatax)/ncase.
            compute rgpmax=cmax(tdatax).
            compute rgpmin=cmin(tdatax).
            compute rgpdes={t(rgpmean),rgpsd,t(rgpmin),t(rgpmax)}.
            compute tdatanm2={nmy,nmx2,'       .'}.
            print rgpdes/title='*************************'+
                '******** Descriptive Statistics ********************************'
                /rname=tdatanm2/clabels='Mean' 'Std.Dev'  'Min' 'Max' /format=!dcms.
            release tdatax.
            @reg3 @rgy=edty/@rgx=tdata1/@rgci=eregci/@rg1vdr=1/@rg1hc=rgprbst
                /@rgALL=1/@rg2cvg=lgconver.
                print /title=!rgline.
                print /title=!rgytitle/format=a8.
            do if (dichy=2).
                print rgycode/title=!rgycode/format=!dcms/rname=prtrslt/cname=rgycodep.
                print rg2ifm/title '>>>>> Model Summary'/format=!dcms/rname=prtrslt/cname=rg2ifmp.
                print rg2class/title '>>>>> '+
    'Classification'/format=!dcms/rname=rg2clspr/cname=rg2clspc.
                print rg2coef(:,1:7)/title '>>>>> '+
    'Coefficients'/format=!dcms/rname=tdatanm/cname=rgbp.
                else.
                print rg1ifm/title '>>>>> Model Summary'/format=!dcms/rname=prtrslt/cname=rg1ifmp.
                compute tdata1={rgbp(1:6),rgbp(:,8:9)}.
                print {rg1coef(:,1:6),rg1coef(:,8:9)}/title '>>>>> '+
    'Coefficients'/format=!dcms/rname=tdatanm/cname=tdata1.
            end if.
            compute rgpemp={0,!rgpemp}.
            compute rgpemq={0,!rgpemq}.
            compute rgpemn1=ncol(rgpemp)-1.
            compute rgpemn2=ncol(rgpemq)-1.
            compute rgpemp=rgpemp(1,2:(rgpemn1+1)).
            compute rgpemq=rgpemq(1,2:(rgpemn2+1)).
            compute rgpemp1=0.
            compute rgpemq1=0.
            compute rgbnvar=nrow(rgb).
            do if (rgbnvar=rgpemn1).
               compute rgpemp1=1.
            end if.
            do if (rgbnvar=rgpemn2).
               compute rgpemq1=1.
            end if.
            do if (rgpemp1=1 or rgpemq1=1).
                compute tdata1=make(1,rgbnvar,0).
                do if (rgpemp1=1).
                    compute tdata1={tdata1;rgpemp}.
                end if.
                do if (rgpemq1=1).
                    compute tdata1={tdata1;rgpemq}.
                    do if (rgpemp1=1).
                        compute tdata1={tdata1;(rgpemp-rgpemq)}.
                    end if.
                end if.
                compute tdata2=nrow(tdata1)-1.
                compute tdata1=tdata1(2:(tdata2+1),:).
                compute rgpempc={'Coeff','S.E.','Z','p','LLCI','ULCI','OR'}. 
                do if(dichy<>2).
                    compute rgpempc(1,3)={'t'}.
                end if.
                do if(rgpemp1=1).
                    compute rgpempr={' Point R'}.
                    do if (rgpemq1=1).
                        compute rgpempr={' Point R';' Point S';' R - S  '}.
                    end if.
                else.
                    compute rgpempr={' Point R'}.
                end if.
                compute tdata3=tdata1*rgb.
                compute tdata3={tdata3,sqrt(diag(tdata1*rgbcov*t(tdata1)))}.
                compute tdata3={tdata3,(tdata3(:,1)&/tdata3(:,2))}.
                do if (dichy=2).
                    compute tdata3={tdata3,2*(1-cdfnorm(abs(tdata3(:,3))))}.
                    compute tdata5=tdata3(:,1).
                    compute tdata5mx=cmax(tdata5).
                    do if (tdata5mx>709.77).
                        compute tdata501=(tdata5>709.77).
                        compute tdata5=tdata5&*(1-tdata501)+tdata501*709.77.
                    end if.
                    compute tdata3={tdata3,(tdata3(:,1)-tdata3(:,2)*ppnd16),
                        (tdata3(:,1)+tdata3(:,2)*ppnd16),exp(tdata5)}.
                    else.
                    compute tdata3={tdata3,2*(1-tcdf(abs(tdata3(:,3)),rg1df2))}.
                    compute tdata3={tdata3,(tdata3(:,1)-tdata3(:,2)*qtout),
                        (tdata3(:,1)+tdata3(:,2)*qtout)}.
                end if.
                print /title '***************************** Marginal Means and Comparison '
                    +'*****************************'.
                print t(tdata1)/title '>>>>> Values plugged into the regression '+
    'model'/rname=tdatanm    
                    /cname=rgpempr/format=!dcms.
                print tdata3/title '>>>>> Point Estimation (and '+
                    'Comparison)'/cname=rgpempc/rname=rgpempr/format=!dcms.
            else.
                print /title=!plotitle.
                print /title=!plonote/space=0.
                !let !ttx1=!head(!e_x)
                !let !ttsta = !concat('Data list free/',!ttx1,!BLANKS(1),!e_y,' SE LLCI ULCI .')
                !let !ttsta = !quote(!ttsta)
                !let !ttend = !concat('Graph/line(multiple)=Value','( ',
                    !e_y,' LLCI ULCI ) by ',!ttx1,'.' )
                !let !ttend =!quote(!ttend)
                compute rgpsd=sqrt(diag(ncase2*(sscp(edtx2)-
                    ((t(csum(edtx2))*csum(edtx2))/ncase)))).
                compute rgpmean=csum(edtx2)/ncase.
                compute rgpxinc=(rgpxmax(1,1)-rgpxmin(1,1))/(rgpdn-1).
                compute rgpdata=make(rgpdn,1,0).
                loop #rgp1=1 to rgpdn.
                    compute rgpdata(#rgp1,1)=rgpxmin(1,1)+rgpxinc*(#rgp1-1).
                end loop.
                do if (dichy=2). 
                    compute qtout=ppnd16.
                end if.
                do if (rgpm>5.5).
                    do if (!rgpmol=3).
                        compute rgpmov={!rgpmov}.
                    else if (!rgpmol=1).
                        compute rgpmov=rgpmean(1,2)-rgpsd(2,1).
                        compute rgpmov={rgpmov,rgpmean(1,2)}.
                        compute rgpmov={rgpmov,rgpmean(1,2)+rgpsd(2,1)}.
                    else if (!rgpmol=2).
                        compute tdata1=edtx2(:,2).
                        @sortv sortv=tdata1/sorttp=1.
                        @prcn prcny=tdata1/prcnp=16.
                        compute rgpmov=percnv.
                        @prcn prcny=tdata1/prcnp=50.
                        compute rgpmov={rgpmov,percnv}.      
                        @prcn prcny=tdata1/prcnp=84.
                        compute rgpmov={rgpmov,percnv}.      
                    end if.
                    compute tdata2=ncol(rgpmov).
                    loop #rgp1=1 to tdata2.
                        print rgpmov(1,#rgp1)/title '====================> '+
                            'Moderator ='/format=!dcms.
                        compute tdata3=make(rgpdn,1,rgpmov(1,#rgp1)).
                        compute tdata1={rgpdata,rgpdata&**2,tdata3,tdata3&**2,
                            rgpdata&*tdata3,rgpdata&**2&*tdata3,rgpdata&*tdata3&**2}.
                        do if (rgpm=6).
                            compute tdata1={tdata1(:,1),tdata1(:,3),tdata1(:,5)}.
                            else if (rgpm=7).
                            compute tdata1={tdata1(:,1:3),tdata1(:,5)}.
                            else if (rgpm=8).
                            compute tdata1={tdata1(:,1:3),tdata1(:,5:6)}.
                            else if (rgpm=9).
                            compute tdata1={tdata1(:,1),tdata1(:,3:5)}.
                            else if (rgpm=10).
                            compute tdata1={tdata1(:,1),tdata1(:,3:5),tdata1(:,7)}.
                        end if.
                        compute tdata1={make(rgpdn,1,1),tdata1}.
                        do if (nvarx2>2).
                            loop #rgp2=3 to nvarx2. 
                                compute tdata1={tdata1,make(rgpdn,1,rgpmean(1,#rgp2))}.
                            end loop.
                        end if.
                        compute rgpdata2={tdata1(:,2),tdata1*rgb,
                            sqrt(diag(tdata1*rgbcov*t(tdata1)))}.
                        compute rgpdata2={rgpdata2,(rgpdata2(:,2)-rgpdata2(:,3)*qtout),
                            (rgpdata2(:,2)+rgpdata2(:,3)*qtout)}.
                        print /title=!ttsta/format=a8/space=1.
                        print rgpdata2/title 'Begin data.'/format=!dcms/space=0.
                        print /title 'End data.'/format=a8/space=0.
                        print /title=!ttend/format=a8/space=0.
                    end loop.
                else.
                    compute tdata1=rgpdata.
                    loop #rgp2=1 to rgpm.
                        compute rgpdata={rgpdata,tdata1&**#rgp2}.
                    end loop.
                    compute rgpdata(:,1)=make(rgpdn,1,1).
                    do if (nvarx2>1).
                        loop #rgp3=2 to nvarx2. 
                            compute rgpdata={rgpdata,make(rgpdn,1,rgpmean(1,#rgp3))}.
                        end loop.
                    end if.
                    compute rgpdata2=make(rgpdn,5,0).
                    compute rgpdata2(:,2)=rgpdata*rgb.
                    compute rgpdata2(:,3)=sqrt(diag(rgpdata*rgbcov*t(rgpdata))).
                    compute rgpdata2(:,1)=rgpdata(:,2).
                    compute rgpdata2(:,4:5)={(rgpdata2(:,2)-rgpdata2(:,3)*qtout),
                        (rgpdata2(:,2)+rgpdata2(:,3)*qtout)}.
                    print /title=!ttsta/format=a8/space=1.
                    print rgpdata2/title 'Begin data.'/format=!dcms/space=0.
                    print /title 'End data.'/format=a8/space=0.
                    print /title=!ttend/format=a8/space=0.
                end if.
            end if.
        end if.
    end if.
end if.
do if (emod=6).
      print {!pcrbtn}/title'Bootstrap'/rname={'Samples:'}/format=f10.0.
   do if (!pcrbtm=1).
      do if (!quote(!seed)="random").
         print {'Random'}/title='Type   : Percentile'/rlabel='Seed   :'/format=a8/space=0.
         else.
         print seed/title='Type   : Percentile'/rlabel='Seed   :'/format=a8/space=0.
      end if.
      else.
      do if (!quote(!seed)="random").
         print {'Random'}/title='Type   : Bias-corrected'/rlabel='Seed   :'/format=a8/space=0.
         else.
         print seed/title='Type   : Bias-corrected'/rlabel='Seed   :'/format=a8/space=0.
      end if.
   end if.
   do if (sero=0).
      do if (nvarx2<2).
         compute ero={ero;5}.
      end if.
       do if (nvarx2>100).
           compute ero={ero;6}.
       end if.
   end if.
    do if (dichy=2).
        compute ero={ero;15}.
    end if.
   compute sero=csum(ero).
   do if (sero=0).
      do if (nvarx2>(ncase-2)).
         compute ero={ero;7}.
      end if.
   end if.
   compute sero=csum(ero).
   do if (sero=0).
      compute kr2b=make(nvarx2,nvarx2+2,-999).
      loop #z= 1 to nvarx2.
         @PCAreg2 pcr_y=edty/pcr_x=edtx2/pcr_n=#z.
         compute kr2b(#z,:)={#z,modelr2,t(stdxcoef)}.
      end loop.
      compute edtplot={kr2b(:,1),evalue,kr2b(:,2:(nvarx2+2))}.
      compute tdata1=t({nmy,nmx2}).
      compute tdata2=t(tdata1).
      do if (!pcrdes=1).
         compute xymimx={cmin({edty,edtx2});cmax({edty,edtx2})}.
         print t({xysdm;xymimx})/title='*************************'+
                '******** Descriptive Statistics ********************************'
               /rname=tdata1/clabels='Mean' 'Std.Dev'  'Min' 'Max' 
    /format=!dcms.
      end if.
        do if (!pcrcor=1).
            print corrx/title='********************************** Correlation Matrix '+
    '**********************************' /cname=tdata2/rname=tdata1/format=!dcms.
        end if.
        do if (!pcrallrb=1).
            compute tdata3={'No.','R-sq',tdata2}.
            print kr2b/title='***************** Number of Components, R-square, and '+
    'Std.Coefficients *****************' /cname=tdata3/rname={'Results:'}/format=!dcms.
        end if.
        compute eva2=evalue/nvarx2*100.
        compute eva3=make(nvarx2,1,0).
        compute evaall={make(nvarx2,1,0),evalue,eva2,eva3}.
        loop #i=1 to nvarx2.
            compute evaall(#i,4)=csum(evalue(1:#i,1))/nvarx2*100.
            compute evaall(#i,1)=#i.
        end loop.
        compute allevanm={'No.','EigenValue','PV','CPV'}.
        do if (!pcreva=1).
            print evaall/title='************************************* Eigenvalues '+
    '**************************************' /cname=allevanm/rname={'Results:'}/format=!dcms.
            print /title 'Notes  :'/format=a8/space=0.
            print /title '1.PV   = Percentage of Variance.'/format=a8/space=0.
            print /title '2.CPV  = Cumulative Percentage of Variance.'/format=a8/space=0.
        end if.
        do if (!pcreve=1).
            compute evaall2={make(nvarx2,1,0),evalue,evector}.
            loop #i=1 to nvarx2.
                compute evaall2(#i,1)=#i.
            end loop.
            compute tdata1={'No.','EigenValue',nmx2}.
            print evaall2/title='************************************** Eigenvectors '+
    '************************************' /cname=tdata1/rname={'Results:'}/format=!dcms.
        end if.
        compute tdata1={'1','2','3','4','5','6','7','8','9','10','11','12','13',
      '14','15','16','17','18','19','20','21','22','23','24','25','26','27','28',
      '29','30','31','32','33','34','35','36','37','38','39','40','41','42','43',
      '44','45','46','47','48','49','50','51','52','53','54','55','56','57','58',
      '59','60','61','62','63','64','65','66','67','68','69','70','71','72','73',
      '74','75','76','77','78','79','80','81','82','83','84','85','86','87','88',
      '89','90','91','92','93','94','95','96','97','98','99','100'}.    
        compute tdata2=tdata1(1,1:nvarx2).
        do if (!pcrcm=1).
            print cpmatrix/title='************************************ Component Matrix '+
    '**********************************' /rname=nmx2/cname=tdata2/format=!dcms.
        end if.
        do if (!pcrcsm=1).
            print scmatrix/title='*************************** Component Score Coefficient'+
            ' Matrix *************************' /rname=nmx2/cname=tdata2/format=!dcms.
        end if.
        do if (!pcrtype=1).
                compute tdata1=1/(ncase-1).
                compute tdata2=make(nvarx2,4,0).
                compute pcrnum=rnd(!pcrnum).        
                do if (pcrnum<1000).
                    compute pcrnum=1000.
                end if.
                compute bteval=make(pcrnum,nvarx2,-999).
                loop #i=1 to pcrnum.
                    compute randmx = sqrt(2*(ln(uniform(ncase,nvarx2))*-1))&*
                        cos(6.28318530717958*uniform(ncase,nvarx2)).
                    compute vcv = tdata1*(sscp(randmx)-((t(csum(randmx))*csum(randmx))/ncase)).
                    compute d = inv(mdiag(sqrt(diag(vcv)))).
                    compute bteval(#i,:)=t(eval(d*vcv*d)).
                end loop.
                loop #i=1 to nvarx2.
                    @sortv sortv=bteval(:,#i)/sorttp=1.
                    @prcn prcny=bteval(:,#i)/prcnp=95.
                    compute tdata2(#i,3)=percnv.
                end loop. 
                compute tdata2(:,2)=evalue.
                loop #i=1 to nvarx2.
                    compute tdata2(#i,1)=#i.
                    do if(tdata2(#i,2)>tdata2(#i,3)).
                        compute tdata2(#i,4)=1.
                    end if. 
                end loop.
                compute pcrpan=csum(tdata2(:,4)).
                do if (pcrpan=0).
                     compute pcrpan=1.
                end if.
                print /title'*********************************** Parallel Analysis '+
                    '**********************************' /format=a8.
                compute pcrnum=!pcrnum.
                do if (pcrnum<1000).
                    !let !tt =!concat('Number of simulations : ','1000')
                    !let !tt=!quote(!tt)
                print tdata2/title=!tt/cname={'No.','Eigenvalue','PA-95th','Extract?'}
                        /rname={'Results:'}/format=!dcms/space=0.
                end if.
                do if (pcrnum>=1000).
                    !let !tt =!concat('Number of simulations : ',!pcrnum)
                    !let !tt=!quote(!tt)
                print tdata2/title=!tt/cname={'No.','Eigenvalue','PA-95th','Extract?'}
                    /rname={'Results:'}/format=!dcms/space=0.
                end if.
                print /title 'Note   : PA-95th = Parallel Analysis-95th '+
                        'percentile.'/format=a8/space=0.
            else if(!pcrtype=2).
                compute pcrnum=!pcrnum.
                do if (pcrnum>100).
                    compute pcrnum=80.
                    else if (pcrnum<50).
                    compute pcrnum=80.
                    else.
                    compute pcrnum=pcrnum.
                end if.
                compute tdata1=0.
                compute pcrpan=0.
                loop #i=1 to nvarx2.
                    do if (evaall(#i,4)>pcrnum and tdata1=0).
                        compute pcrpan=#i.
                        compute tdata1=1.
                    end if.
                    do if (#i=nvarx2).
                        do if (pcrpan=0).
                            compute pcrpan=nvarx2.
                        end if.
                    end if.
                end loop.
            else if(!pcrtype=3).
                compute pcrnum=!pcrnum.
                do if (pcrnum>nvarx2).
                    compute pcrnum=1.
                    else if (pcrnum<0).
                    compute pcrnum=1.
                    else.
                    compute pcrnum=pcrnum.
                end if.
                compute pcrpan=0.
                loop #i=1 to nvarx2.
                    do if(evaall(#i,2)>pcrnum).
                        compute pcrpan=pcrpan+1.
                    end if.
                end loop.
            else if(!pcrtype=4).
                compute pcrnum=rnd(!pcrnum).
                do if (pcrnum>nvarx2).
                    compute pcrnum=1.
                    else if (pcrnum<1).
                    compute pcrnum=1.
                end if.
                compute pcrpan=pcrnum.
        end if.
        compute pcrbtn=!pcrbtn.
        compute pcrcoef1=make(pcrbtn,(nvarx2+1),-999).
        compute pcrcoef2=make(pcrbtn,nvarx2,-999).
        compute btnero=1.
        compute pcrbtdt={edty,edtx2}.
        compute rawdata=pcrbtdt.
        compute tdataq=0.
        loop #btn=1 to pcrbtn.
            loop #btsp=1 to ncase.
                compute pcrbtdt(#btsp,:)=rawdata((trunc(uniform(1,1)*ncase)+1),:).
            end loop.            
            @candr crdata=pcrbtdt.
            do if (candrn>0).
                compute btero=btero+1.
                else.
                @PCAreg2 pcr_y=pcrbtdt(1:ncase,1)/pcr_x=pcrbtdt(1:ncase,2:(nvarx2+1))
                    /pcr_n=pcrpan.
                compute pcrcoef1(btnero,:)=t(xcoef).
                compute pcrcoef2(btnero,:)=t(stdxcoef).
                compute btnero=btnero+1.
            end if.
        end loop.
        do if (btero>0).
            loop #i=1 to btero.
                compute tdata1=trunc(uniform(1,1)*btero)+1.
                compute pcrcoef1((pcrbtn+1-#i),:)=pcrcoef1(tdata1,:).
                compute pcrcoef2((pcrbtn+1-#i),:)=pcrcoef2(tdata1,:).
            end loop.
        end if.
        @PCAreg2 pcr_y=edty/pcr_x=edtx2/pcr_n=pcrpan.
        loop #i=1 to (nvarx2+1).
            @sortv sortv=pcrcoef1(:,#i).
            do if (#i>1).
                @sortv sortv=pcrcoef2(:,(#i-1)).
            end if.
        end loop.
        @btsot btrawd=t(xcoef)/btalld=pcrcoef1/btcitp=pcrbtm
            /btnewnm=xdm1/btcih=eregcih/btcil=eregcil.
        @btsot btrawd=t(stdxcoef)/btalld=pcrcoef2/btcitp=pcrbtm
            /btnewnm=xdm2/btcih=eregcih/btcil=eregcil.
        print /title'***************************** Principal Component Rgression '+
    '****************************'/format=a8.
        print pcrpan/title 'Number of components:'/format=f8.0/space=0.
        compute pcrsmry={sqrt(modelr2),modelr2}.
        compute pcrsmry={pcrsmry,(1-(1-pcrsmry(1,2))*(ncase-1)/(ncase-1-nvarx2))}.
        print pcrsmry/title '>>>>> Model Summary'/rname={'Results:'}
                /cname={'R','R-square','Adj-Rsq'}/format=!dcms. 
        print ftest/title '>>>>> ANOVA' /cname={'SS','df','MS','F','p'}
                /rname={'Model','Error','Total','       .'}/format=!dcms.
        compute xdm1=t(xdm1).
        compute xdm2=t(xdm2).
        compute tdata1=xdm1(:,1)&/xdm1(:,4).
        compute tdata2=2*(1-tcdf(abs(tdata1),ftest(2,2))).
        compute regnm={'Constant';t(nmx2)}.
        compute tdata3={0;xdm2(:,1)}.
        print {xdm1(:,1:2),xdm1(:,4),tdata1,tdata2,xdm1(:,5:6),tdata3}
            /title '>>>>> Coefficients'
            /cname={'Coeff','BootMean','BootSE','t','p','BootLLCI','BootULCI','Std.Coef'}
            /rname=regnm/format=!dcms.
        compute tdata1=xdm2(:,1)&/xdm2(:,4).
        compute tdata2=2*(1-tcdf(abs(tdata1),ftest(2,2))).
        compute tdata1=!pcrplot1+!pcrplot2+!pcrplot3.
        do if (tdata1>0).
            print /title=!plotitle.
            print /title=!plonote/space=0.
            do if (!pcrplot3=1).
                print /title'====================> Number of components vs. Eigenvalues'.
                print /title 'Data list free/Component Eigenvalue.'/format=a8.
                print {edtplot(:,1),edtplot(:,2)}/title 'Begin data.'/format=!dcms/space=0.
                print /title 'End data.'/format=a8/space=0.
                print /title 'Graph/scatterplot=Component with Eigenvalue .'/foramt=a8/space=0.
            end if.
            do if (!pcrplot1=1).
                print /title'====================> Number of components vs. R-squares'.
                print /title 'Data list free/Component Rsquare.'/format=a8.
                print {edtplot(:,1),edtplot(:,3)}/title 'Begin data.'/format=!dcms/space=0.
                print /title 'End data.'/format=a8/space=0.
                print /title 'Graph/scatterplot=Component with Rsquare .'/foramt=a8/space=0.
            end if.
            do if (!pcrplot2=1).
                print /title'====================> Number of components vs. Std.Coefficients'.
                !let !tt = !concat('   Component')
                !do !vars !in (!e_x)
                    !let !tt = !concat(!tt,!blanks(1),!vars)
                !doend
                !let !tt = !concat('Data list free/',!tt)
                !let !tt = !concat(!tt,' .')
                !let !tt =!quote(!tt)
                print /title=!tt /foramt=a8.
                print {edtplot(:,1),edtplot(:,4:(nvarx2+3))}
                    /title 'Begin data.'/format=!dcms/space=0.
                print /title 'End data.'/format=a8/space=0.
                !let !tt = !concat('Graph/line(multiple)=Value','(')
                !do !vars !in (!e_x)
                    !let !tt = !concat(!tt,!blanks(1),!vars)
                !doend
                !let !tt = !concat(!tt,' ) BY Component.')
                !let !tt =!quote(!tt)
                print /title=!tt/format=a8/space=0.
            end if.
        end if.
   end if.
end if.
do if (emod=8).
   compute spltp=rnd(!spltp).
   compute splktp=rnd(!splktp).
   compute splqv={0,!splqv}.
   compute splrbst=!splrbst.
   compute splsave=rnd(!splsave).
   do if (dichy<>2).
       do if (splrbst<>9).
           !let !tt = !concat('Standard errors for regression coefficients : HC',!splrbst)
           !let !tt = !quote (!tt)
               print /title=!tt/format=a8.
       end if .
   end if.
   do if (spltp<4).
      compute splkvm=make(99,99,0).
      loop #i=1 to 99.
         compute tdata1=#i+1.
         compute tdata2=100/tdata1.
         loop #i2=2 to tdata1.
            compute splkvm(#i,(#i2-1))=(#i2-1)*tdata2.
         end loop.
      end loop.
      else.
      compute splkvm=make(99,99,0).
      compute splkvm(3:6,1:7)={10,50,90,0,0,0,0;5,35,65,95,0,0,0;
          5,27.5,50,72.5,95,0,0;5,23,41,59,77,95,0}.
      loop #rcs=7 to 99.
          compute tdata1=95/(#rcs-1).
          loop #rcs2=1 to #rcs.
              compute splkvm(#rcs,#rcs2)=2.5+(#rcs2-1)*tdata1.
          end loop.
      end loop.
   end if.
   compute splxnm={"Spline1","Spline2","Spline3","Spline4","Spline5","Spline6","Spline7",
       "Spline8","Spline9","Spline10","Spline11","Spline12","Spline13","Spline14","Spline15",
       "Spline16","Spline17","Spline18","Spline19","Spline20","Spline21","Spline22","Spline23",
       "Spline24","Spline25","Spline26","Spline27","Spline28","Spline29","Spline30","Spline311",
       "Spline32","Spline33","Spline34","Spline35","Spline36","Spline37","Spline38","Spline39",
       "Spline40","Spline41","Spline42","Spline43","Spline44","Spline45","Spline46","Spline47",
       "Spline48","Spline49","Spline50","Spline51","Spline52","Spline53","Spline54","Spline55",
       "Spline56","Spline57","Spline58","Spline59","Spline60","Spline61","Spline62","Spline63",
       "Spline64","Spline65","Spline66","Spline67","Spline68","Spline69","Spline70","Spline71",
       "Spline72","Spline73","Spline74","Spline75","Spline76","Spline77","Spline78","Spline79",
       "Spline80","Spline81","Spline82","Spline83","Spline84","Spline85","Spline86","Spline87",
       "Spline88","Spline89","Spline90","Spline91","Spline92","Spline93","Spline94","Spline95",
       "Spline96","Spline97","Spline98","Spline99"}.
   compute splxnm2=splxnm.
   compute splkn=ncol(splqv).
   do if (splkn=1).
      compute ero={ero;9}.
   end if.
   compute sero=csum(ero).
   do if (sero=0).
      compute splqv=splqv(1,2:splkn).
      do if (splktp=1).
          compute splkn=splqv(1,1).
          else.
          compute splkn=splkn-1.
      end if.
      do if (spltp<4).
         do if (splktp=1).
            compute tdata1=rnd(splqv(1,1)).
            do if (tdata1<1 or tdata1>99).
               compute ero={ero;10}.
            end if.
            compute sero=csum(ero).
            do if (sero=0).
               compute splqv2=t(splkvm(tdata1,1:tdata1)).
               compute tdata2=edtx1.
               @sortv sortv =tdata2/sorttp=1.
               compute splkvm2=make(tdata1,1,0).
               loop #i = 1 to tdata1.
                  @prcn prcny=tdata2/prcnp=splqv2(#i,1).
                  compute splkvm2(#i,1)=percnv.
               end loop.
            end if.
         else.
            do if (splkn<1 or splkn>99).
               compute ero={ero;10}.
            end if.
            compute tdata1=rmax(splqv).
            compute tdata2=rmin(splqv).
            do if (tdata1>=edtx1max or tdata2<=edtx1min).
               compute ero={ero;11}.
            end if.
            compute sero=csum(ero).
            do if (sero=0).
               compute splkvm2=t(splqv).
               @sortv sortv=splkvm2/sorttp=1.
            end if.
         end if.
      else.
         do if (splktp=1).
            compute tdata1=rnd(splqv(1,1)).
            do if (tdata1<3 or tdata1>99).
               compute ero={ero;13}.
            end if.
            compute sero=csum(ero).
            do if (sero=0).
               compute splqv2=t(splkvm(tdata1,1:tdata1)).
               compute tdata2=edtx1.
               @sortv sortv =tdata2/sorttp=1.
               compute splkvm2=make(tdata1,1,0).
               loop #i = 1 to tdata1.
                  @prcn prcny=tdata2/prcnp=splqv2(#i,1).
                  compute splkvm2(#i,1)=percnv.
               end loop.
            end if.
         else.
            do if (splkn<3  or splkn>99).
               compute ero={ero;13}.
            end if.
            compute tdata1=rmax(splqv).
            compute tdata2=rmin(splqv).
            do if (tdata1>=edtx1max or tdata2<=edtx1min).
               compute ero={ero;11}.
            end if.
            compute sero=csum(ero).
            do if (sero=0).
               compute splkvm2=t(splqv).
               @sortv sortv=splkvm2/sorttp=1.
            end if.
         end if.
      end if. 
      compute sero=csum(ero).
      do if (sero=0).
         compute splplot=!splplot.
         compute splplotn=!splplotn.
         compute splemmp={!splemmp}.
         compute splemmn=ncol(splemmp)-1.
         compute splemtyp=0.
         compute splmean=csum(edtx2)/ncase.
         compute splmax=cmax(edtx2).
         compute splmin=cmin(edtx2).
         do if (splemmn>0).
             compute splemmp=splemmp(1,2:(splemmn+1)).
             compute tdata1=0.
             loop #sple=1 to splemmn.
                 do if (splemmp(1,#sple)>=splmin(1,1) and splemmp(1,#sple)<=splmax(1,1)).
                     compute tdata1={tdata1;splemmp(1,#sple)}.
                 end if.
             end loop.
             compute splplotn=nrow(tdata1)-1.
             do if (splplotn>0).
                compute splemmp=tdata1(2:(splplotn+1),1).
                compute splemtyp=1.
             end if.
             do if (splplotn=2).
                compute splemtyp=2.
             end if.
         else .
             do if (splplot=1).
                compute splemtyp=3.
                compute splinc=(splmax(1,1)-splmin(1,1))/(splplotn-1).
                loop #splp=1 to splplotn.
                    do if (#splp=1).
                        compute splemmp=splmin(1,1).
                    else.
                        compute splemmp={splemmp;(splmin(1,1)+splinc*(#splp-1))}.
                    end if.
                end loop.
             end if.
         end if.
         do if (splktp=1).
            do if (spltp<4).
                compute tdata1=splxnm(1,1:splkn).
                print {t(splqv2);t(splkvm2)}/title '>>>>> Percentiles and Knot Values'
                      /rname={'PCTL','Value','       .'}/cname=tdata1/format=!dcms.
                else.
                print {t(splqv2);t(splkvm2)}/title '>>>>> Percentiles and Knot Values'
                      /rname={'PCTL','Value','       .'}/format=!dcms.
            end if.
            else.
            do if (spltp<4).
                compute tdata1=splxnm(1,1:splkn).
                print t(splkvm2)/title 'Knot Values'/rname={'Value','       .'}
                /cname=tdata1/format=!dcms.
            else.
                print t(splkvm2)/title 'Knot Values'/rname={'Value','       .'}/format=!dcms.
            end if.
         end if.
         do if (splemtyp>0).
             compute edtx1={edtx1;splemmp}.
             compute ncase=ncase+splplotn.
         end if.
         do if (spltp<4).
            compute splx2=make(ncase,splkn,0).
            compute splx2s=make(ncase,splkn,0).
            compute tdata1=make(ncase,1,spltp).
            loop #i= 1 to splkn.
                compute splx2(:,#i)=edtx1-splkvm2(#i,1).
                compute splx2s(:,#i)=(edtx1>splkvm2(#i,1)).
                compute splx2(:,#i)=(splx2(:,#i)&*splx2s(:,#i))&**spltp.
            end loop.
            else .
            compute splx2=make(ncase,splkn-2,0).
            compute splkend=splkvm2(splkn,1).
            compute splkend2=splkvm2((splkn-1),1).
        loop #i=1 to (splkn-2).
            compute tdata1=(edtx1-splkvm2(#i,1)).
            compute tdata2=(edtx1-splkend2).
            compute tdata3=(edtx1-splkend).
            compute tdata1=(tdata1&*(tdata1>0))&**3.
            compute tdata2=(tdata2&*(tdata2>0))&**3*(splkend-splkvm2(#i,1))/(splkend-splkend2).
            compute tdata3=(tdata3&*(tdata3>0))&**3*(splkend2-splkvm2(#i,1))/(splkend-splkend2). 
            compute splx2(:,#i)=tdata1-tdata2+tdata3.
        end loop.
            compute splrcsc=1.
            do if (splrcsc=1).
                compute splx2=splx2/((splkend-splkvm2(1,1))**2).
            end if.
         end if.
         do if (splemtyp>0).
             compute edtx1=edtx1(1:(ncase-splplotn),:).
             compute splemx2=splx2((ncase-splplotn+1):ncase,:).
             compute splx2=splx2(1:(ncase-splplotn),:).
             compute ncase=ncase-splplotn.
             compute splemmp={make(splplotn,1,1),splemmp}.
         end if.
         do if (spltp=1).
             compute splxall={edtx1,splx2}.
             comptue splxnm={nmx1,splxnm(1,1:splkn)}.
             do if (splemtyp>0).
                  compute splxallp={splemmp,splemx2}.
             end if.
         else if (spltp=2).
             compute splxall={edtx1,edtx1&**2,splx2}.
             comptue splxnm={nmx1,'X^2',splxnm(1,1:splkn)}.
             do if (splemtyp>0).
                  compute splxallp={splemmp,splemmp(:,2)&**2,splemx2}.
             end if.
         else if (spltp=3).
             compute splxall={edtx1,edtx1&**2,edtx1&**3,splx2}.
             comptue splxnm={nmx1,'X^2','X^3',splxnm(1,1:splkn)}.
             do if (splemtyp>0).
                  compute splxallp={splemmp,splemmp(:,2)&**2,splemmp(:,2)&**3,splemx2}.
             end if.
         else if (spltp=4).
             compute splxall={edtx1,splx2}.
             comptue splxnm={nmx1,splxnm(1,1:(splkn-2))}.
             do if (splemtyp>0).
                  compute splxallp={splemmp,splemx2}.
             end if.
         end if.
         do if (nvarx2>1).
            compute splxall={splxall,edtx2(:,2:nvarx2)}.
            compute splxnm={splxnm,nmx2(1,2:nvarx2)}.
            do if (splemtyp>0).
                loop #spl2=2 to nvarx2.
                    compute splxallp={splxallp,make(splplotn,1,splmean(1,#spl2))}.
                end loop.
            end if.
         end if.
         compute splxnm={'Constant',splxnm,'       .'}.
         @candr crdata=splxall.
         compute ero={ero;candrn}.
         compute sero=csum(ero).
         do if (sero=0).
            print /title=!rgline.
            print /title=!rgytitle/format=a8.
            @reg3 @rgy=edty/@rgx=splxall/@rgci=eregci/@rg1vdr=1/@rg1hc=splrbst
                /@rgALL=1/@rg2cvg=lgconver.
            compute splicnm={'AIC','BIC','AICc','aBIC'}.
            do if (dichy=2).
                print rgycode/title=!rgycode/format=!dcms/rname=prtrslt/cname=rgycodep.
                print rg2ifm/title '>>>>> Model Summary'/format=!dcms
                     /rname=prtrslt/cname=rg2ifmp.
                print rgics/title '>>>>> Information Criterion'/rname=prtrslt
                    /cname=rgicsp/format=!dcms.
                print rg2class/title '>>>>> Classification'/format=!dcms
                    /rname=rg2clspr/cname=rg2clspc.
                print rg2coef(:,1:7)/title '>>>>> Coefficients'/format=!dcms
                    /rname=splxnm/cname=rgbp.
                else.
                print rg1ifm/title '>>>>> Model Summary'/format=!dcms
                    /rname=prtrslt/cname=rg1ifmp.
                compute aic=ncase*ln(rg1ess/ncase)+2*(rgnvar+1).
                compute bic=ncase*ln(rg1ess/ncase)+(rgnvar+1)*ln(ncase).
                compute aicc=aic+2*(rgnvar+1)*(rgnvar+2)/(ncase-(rgnvar+1)-1).
                compute abic=ncase*ln(rg1ess/ncase)+ln((ncase+2)/24)*(rgnvar+1).
                print rgics/title '>>>>> Information Criterion'/rname=prtrslt
                    /cname=rgicsp/format=!dcms.
                compute tdata1={rgbp(:,1:6),rgbp(:,8:9)}.
                print {rg1coef(:,1:6),rg1coef(:,8:9)}/title '>>>>> Coefficients'/format=!dcms
                    /rname=splxnm/cname=tdata1.
            end if.
            do if (splemtyp>0).
                print /title '************************************* Marginal Means '+
                    '************************************'.
                print /title 'Covariates appearing in the model are evaluated at mean '+
                    'value.'/format=a8/space=0.
                compute splempc={nmx1,'Coeff','S.E.','t','p','LLCI','ULCI','OR'}. 
                do if(dichy=2).
                    compute splempc(1,3)={'Z'}.
                    compute qtout=ppnd16.
                end if.
                do if (splemtyp=2).
                    compute splxallp={splxallp;splxallp(1,:)-splxallp(2,:)}.
                    compute splemmp={splemmp;splemmp(1,:)-splemmp(2,:)}.
                end if.
                compute tdata1=splxallp*rgb.
                compute tdata1={tdata1,sqrt(diag(splxallp*rgbcov*t(splxallp)))}.
                compute tdata1={tdata1,(tdata1(:,1)&/tdata1(:,2))}.
                do if (dichy=2).
                    compute tdata1={tdata1,2*(1-cdfnorm(abs(tdata1(:,3))))}.
                    else.
                    compute tdata1={tdata1,2*(1-tcdf(abs(tdata1(:,3)),rg1df2))}.
                end if.
                compute tdata1={splemmp(:,2),tdata1,tdata1(:,1)-tdata1(:,2)*qtout,
                    tdata1(:,1)+tdata1(:,2)*qtout}.
                do if (splemtyp=1).
                    print tdata1/title'====================> Point Estimation'/rname=prtrslt
                        /cname=splempc/format=!dcms.
                else if (splemtyp=2).
                    compute tdata2={'Point R:','Point S:',' R - S :'}.
                    print tdata1/title '====================> Point Estimation (and Comparison)'
                        /rname=tdata2/cname=splempc/format=!dcms.
                else if (splemtyp=3).
                    print /title '====================> Code for Plot'/format=!dcms.
                    print /title=!plonote/space=0.
                    !let !ttx1=!head(!e_x)
                    !let !ttsta = !concat('Data list free/',!ttx1,!BLANKS(1),!e_y,' SE LLCI ULCI '+
    '.')    
                    !let !ttsta = !quote(!ttsta)
                    !let !ttend = !concat('Graph/line(multiple)=Value','( ',
                        !e_y,' LLCI ULCI ) by ',!ttx1,'.' )
                    !let !ttend =!quote(!ttend)
                    print /title=!ttsta/format=a8/space=1.
                    print {tdata1(:,1:3),tdata1(:,6:7)}/title 'Begin data.'/format=!dcms/space=0.
                    print /title 'End data.'/format=a8/space=0.
                    print /title=!ttend/format=a8/space=0.
                end if.
            end if.
         end if.
        do if (splsave=1).
            print /title '*************************************** Save Data '+
            '***************************************'.
            do if (spltp=1).
                compute savedt={eregid,edty,edtx2,splx2}.
                compute splkn=splkn.
                else if (spltp=2).
                compute savedt={eregid,edty,edtx2,edtx1&**2,splx2}.
                compute splxnm2={'X^2',splxnm2}.
                compute splkn=splkn+1.
                else if (spltp=3).
                compute savedt={eregid,edty,edtx2,edtx1&**2,edtx1&**3,splx2}.
                compute splxnm2={'X^2','X^3',splxnm2}.
                compute splkn=splkn+2.
                else if (spltp=4).
                compute savedt={eregid,edty,edtx2,splx2}.
                compute splkn=splkn-2.
            end if.
            compute nvarsv=nvarx2+2.
            save savedt/outfile=*/vars=EReg_ID !e_y !e_x.
            compute prnm={colnm(1,(nvarsv+1):(nvarsv+splkn));splxnm2(1,1:splkn)}. 
            print t(prnm)/title'>>>>> Column Names vs. Variables' 
                /cname={'COL Name','Variable'}/rname={'List   :'}/format=a8.
        end if.
      end if.
   end if.
end if.
do if (emod=2).
    compute mrdcasen=!mrdcasen.
    compute mrassot=!mrassot.
    do if (dichy<>2).
        compute mrrbste=!mrrbste.
        do if (mrrbste<>9).
            !let !tt = !concat('Standard errors for regression coefficients : HC',!mrrbste)
            !let !tt = !quote (!tt)
                print /title=!tt/format=a8.
        end if .
    end if.
    do if (nvarx2>(ncase-2)).
         compute ero={ero;7}.
    end if.
   compute sero=csum(ero).
    do if (sero=0).
        compute tdata2={edty,edtx2}.
        compute mrvcv=ncase2*(sscp(tdata2)-((t(csum(tdata2))*csum(tdata2))/ncase)).
        compute mrd=inv(mdiag(sqrt(diag(mrvcv)))).
        compute mrcr=mrd*mrvcv*mrd.
        compute mrmean=t(csum(tdata2)/ncase).
        compute mrsd=sqrt(diag(mrvcv)).
        do if (!mrdes=1).
            compute tdata1={mrmean,mrsd,t(cmin(tdata2)),t(cmax(tdata2))}.
            compute tdata2=t({nmy,nmx2}).
            print tdata1/title='********************************* Descriptive '
                +'Statistics ********************************'
               /rname=tdata2 /clabels='Mean' 'Std.Dev'  'Min' 'Max' /format=!dcms.
        end if.
        print /title=!rgline.
        print /title=!rgytitle.
       @reg3 @rgy=edty/@rgx=edtx2/@rgci=eregci/@rg1vdr=1/@rg1hc=mrrbste
                /@rgALL=1/@rg2cvg=lgconver.
        compute tdatanm={'Constant',nmx2}.
        do if (dichy=2).
            compute rg2hstrp={'Iteration','-2LL',tdatanm}.
            print rgycode/title=!rgycode/format=!dcms/rname=prtrslt/cname=rgycodep.
            print rg2hstr/title '>>>>> Iteration History'/format=!dcms
                /rname=prtrslt/cname=rg2hstrp.
            print rg2ifm/title '>>>>> Model Summary'/format=!dcms
                /rname=prtrslt/cname=rg2ifmp.
            print rgics/title '>>>>> Information Criterion'/rname=prtrslt
                /cname=rgicsp/format=!dcms.
            print rg2class/title '>>>>> Classification'/format=!dcms
                /rname=rg2clspr/cname=rg2clspc.
            print rg2coef/title '>>>>> Coefficients'/format=!dcms
                /rname=tdatanm/cname=rgbp.
        else.
            compute cacudata=rg1coef(2:(nvarx2+1),9).
            compute mrrsq=rg1rsq.
            compute tdata7=rg1predy.
            compute tdata8=rg1e.
            compute tdata9=rg1e&**2.
            compute tdata1=0.
            loop #q1=1 to (ncase-1).
                compute tdata1=tdata1+(tdata8((#q1+1),1)-tdata8(#q1,1))**2.
            end loop.
            compute rg1smry1={rg1ifm(1,1)**0.5,rg1ifm(1,1:2),
                rg1ifm(1,3)**0.5,(tdata1/csum(tdata9))}.
            compute rg1smry2={rg1rss,nvarx2,rg1rssm,rg1f,rg1pm;
                                          rg1ess,rg1df2,rg1essm,0,0;
                                          rg1yss,(ncase-1),0,0,0}.
            print rg1smry1/title '>>>>> Model Summary'/format=!dcms/rname=prtrslt
                    /cname={'R','R-sq','Adj.R-sq','RMSE','D-W'}.
            print rgics/title '>>>>> Information Criterion'/rname=prtrslt
                /cname=rgicsp/format=!dcms.
            print rg1smry2/title '>>>>> ANOVA' /foramt=!dcms
                /rname={'Model';'Error';'Total';'       .'}/cname={'SS','df','MS','F','p'}.
            print rg1coef/title '>>>>> Coefficients'/format=!dcms/rname=tdatanm/cname=rgbp.
        end if.
        do if (!mrcovb=1). /** Coefficient matrix.
            compute rgbd=inv(mdiag(sqrt(diag(rgbcov)))).
            compute rgbcr=rgbd*rgbcov*rgbd.
            print /title'********************************** Coefficient Matrix '+
                '***********************************'.
            print rgbcov/title '>>>>> Covariance Matrix'/cname=tdatanm
                /rname=tdatanm/format=!dcms.
            print rgbcr/title '>>>>> Correlation Matrix'/cname=tdatanm
                /rname=tdatanm/format=!dcms.
        end if.
        do if (dichy<>2 and !mrcoefev=1). /** Coefficient evaluaion for OLS.
            compute partr=rg1coef(:,9)&**0.5.
            compute partialr=sqrt(rg1t&**2&/(rg1t&**2+rg1df2))&*(rgb&/abs(rgb)).
            compute coeffeva={mrcr(2:(nvarx2+1),1),partialr(2:(nvarx2+1),1),
                            partr(2:(nvarx2+1),1),rg1coef(2:(nvarx2+1),7)}.
            compute coefevnm={'r','Prtl.r','Semi.r','Std.Coef','Str.Coef','Str.C-sq','RIW'}.
            compute tdata3={coeffeva(:,1)/rg1ifm(1,1)**0.5}.
            compute coeffeva={coeffeva,tdata3,tdata3&**2}.
            compute mrrxx=mrcr(2:(nvarx2+1),2:(nvarx2+1)).
            compute mrrxy=mrcr(2:(nvarx2+1),1).
            call eigen(mrrxx,mrevc,mreig).
            do if (nvarx2>1).
                compute mrlbd1=mrevc*sqrt(mdiag(mreig))*t(mrevc).
                compute tdata2=mrlbd1&**2*(inv(mrlbd1)*mrrxy)&**2.
                compute coeffeva={coeffeva,tdata2}.
            end if.
            print coeffeva/title '******************************** Coefficients '+
                'Evaluation ********************************'/cname=coefevnm
                /rname=nmx2//format=!dcms.
            print /title'Notes   :'/foramt=a8/space=0.
            print /title'1.r        = Pearson Correlation'/format=a8/space=0.
            print /title'2.Prtl.r   = Partial Correlation'/format=a8/space=0.
            print /title'3.Semi.r   = Semipartial Correlation (Part Correlation)'
                /format=a8/space=0.    
            print /title'4.Std.Coef = Standardized Regression Coefficient'/format=a8/space=0.
            print /title'5.Str.Coef = Structure Coefficient'/format=a8/space=0.
            print /title'6.Str.C-sq = Square of Structure Coefficient'/format=a8/space=0.
            do if (nvarx2>1).
               print /title'7.RIW      = Relative Importance Weight'/format=a8/space=0.
            end if.
        end if.
        compute tdata1=!mrres.
        compute tdata2=!mrdsts.
        compute tdata0=tdata1+tdata2.
        do if (dichy<>2 and tdata0>0). /** distance and residuals for OLS.
            print /title '********************************** Diagnosis for Cases '+
                    '**********************************'/format=a8.
            compute rgx={make(ncase,1,1),edtx2}.
            /** distance **/.
            compute rg1hat2=diag(rgx*inv(t(rgx)*rgx)*t(rgx)).
            compute rg1cook=rg1e&**2&*rg1hat2&/((1-rg1hat2)&**2*(nvarx2+1)*rg1essm).
            compute rg1maha=(rg1hat2-1/ncase)*(ncase-1).
            /** residuals **/.
            compute rg1tres=rg1e/sqrt(rg1essm*(1-rg1hat2)).
            compute tdata1=make(ncase,1,(ncase-nvarx2-2)).
            compute tdata2=tdata1+1.
            compute rg1tres2=((tdata1&/(tdata2-rg1tres&**2))&**0.5)&*rg1tres.
            compute stdres=rg1e/rg1essm**0.5.
            /** further details.
            compute stdresp=(1-cdfnorm(abs(stdres)))*2.
            compute rg1tresp=2*(1-tcdf(abs(rg1tres2),rg1df2)).
            compute rg1mahap=1-chicdf(rg1maha,nvarx2).
            compute rg1mhf=ncase*rg1maha*rg1df2
                /(nvarx2*((ncase-1)**2-rg1maha*ncase)).
            compute rg1mhfp=1-fcdf(rg1mhf,nvarx2,rg1df2).
            compute ncaseout=rnd(!mrdcasen).
            do if (ncaseout<5).
                compute ncaseout=5.
            end if.
            do if (ncaseout>ncase).
                compute ncaseout=5.
            end if.
            do if (!mrres=1).
                @sortv sortv=stdresp/sorttp=1.
                compute stdres(sort2)=stdres.
                compute stdresn=eregid.
                compute stdresn(sort2)=eregid.
                @sortv sortv=rg1tresp/sorttp=1.
                compute rg1tres2(sort2)=rg1tres2.
                compute rg1tresn=eregid.
                compute rg1tresn(sort2)=eregid.
                compute tdata1={'EReg_ID','Std.Res','p'}.
                print {stdresn(1:ncaseout,1),stdres(1:ncaseout,1),stdresp(1:ncaseout,1)}
                    /title'====================> Standardized Residuals (Std.Res)'
                    /rname={'Results:'}/cname=tdata1/format=!dcms.
                compute tdata1={'EReg_ID','Stud.Res','df','p'}.
                print {rg1tresn(1:ncaseout,1),rg1tres2(1:ncaseout,1),make(ncaseout,1,rg1df2),
                    rg1tresp(1:ncaseout,1)}
                    /title'====================> Studentized Residuals deleted (Stud.Res)'
                    /rname={'Results:'}/cname=tdata1/format=!dcms.
            end if.
            do if (!mrdsts=1).
                @sortv sortv=rg1maha/sorttp=-1.
                compute rg1mahan=eregid.
                compute rg1mahan(sort2)=eregid.
                compute rg1mhp=rg1mahap.
                compute rg1mhp(sort2)=rg1mahap.
                compute rg1mhf2=rg1mhf.
                compute rg1mhf2(sort2)=rg1mhf.
                compute rg1mhfp2=rg1mhfp.
                compute rg1mhfp2(sort2)=rg1mhfp.
                compute tdata2={rg1mahan,rg1maha,make(ncase,1,nvarx2),rg1mhp,
                rg1mhf2,make(ncase,1,nvarx2),make(ncase,1,rg1df2),rg1mhfp2}.
                compute tdata2=tdata2(1:ncaseout,:).
                compute tdata1={'EReg_ID','MD','Chi.df','Chi.p','F.value','F.df1','F.df2','F.p'}.
                print tdata2 /rname={'Results:'}/cname=tdata1/format=!dcms
                    /title'====================> squared Mahalanobis Distance (MD)'.
                @sortv sortv=rg1cook/sorttp=-1.
                compute rg1ckn=eregid.
                compute rg1ckn(sort2)=eregid.
                compute tdata1={'EReg_ID','CD'}.
                print {rg1ckn(1:ncaseout,1),rg1cook(1:ncaseout,1)}
                    /title"====================> Cook's Distance (CD)"
                    /rname={'Results:'}/cname=tdata1/format=!dcms.
                @sortv sortv=rg1hat2/sorttp=-1.
                compute rg1ckn=eregid.
                compute rg1ckn(sort2)=eregid.
                compute tdata1={'EReg_ID','Leverage','2p/n','3p/n'}.
                compute tdata2=make(ncaseout,1,((nvarx2+1)/ncase*2)).
                compute tdata3=make(ncaseout,1,((nvarx2+1)/ncase*3)).
                print {rg1ckn(1:ncaseout,1),rg1hat2(1:ncaseout,1),tdata2,tdata3}
                    /title"====================> Leverage"
                    /rname={'Results:'}/cname=tdata1/format=!dcms.
                print /title 'Notes  : '/space=0.
                print /title '1.p    = number of IVs + 1.'/space=0.
                print /title '2.n    = sample size.'/space=0.
                print /title '3.SPSS reports the centered leverage value '+
                    'that is equal to leverage - 1/n.'/space=0.
            end if.
        end if.
        do if (dichy<>2 and !mrheter=1). /** Heteroskedasticity Tests.
            print /title'******************************** Heteroskedasticity Tests '+
                    '*******************************'/foramt=a8.
            compute tdata1=nvarx2*2+(nvarx2-1)*nvarx2/2+2.
            do if (ncase<tdata1).
                print /title'>>>>> White Test'/format=a8.
                print/title 'Note   : This test could not be performed due to '+
                    'the low sample size.'/format=a8/space=0. 
            else.
                compute tdata1=edtx2.
                loop #i = 1 to nvarx2.
                    compute tdata1(:,#i)=(edtx2(:,#i)-mrmean((#i+1),1))
                        /sqrt(mrvcv((#i+1),(#i+1))).
                end loop.
                loop #i= 1 to nvarx2.
                    compute tdata1={tdata1,tdata1(:,#i)&*tdata1(:,#i)}.
                    loop #i2=1 to nvarx2.
                        do if (#i2>#i).
                            compute  tdata1={tdata1,tdata1(:,#i)&*tdata1(:,#i2)}.
                        end if.
                    end loop.
                end loop.
                @candr crdata=tdata1.
                do if (candrn=0).
                    @reg3 @rgy=tdata9/@rgx=tdata1/@rgci=eregci/@rg1vdr=0
                                /@rg1hc=9/@rgALL=1/@rg2cvg=lgconver.
                    compute mrht1s=ncase*rg1rsq.
                    compute tdata3={mrht1s,rg1df1,(1-chicdf(mrht1s,rg1df1))}.
                    print tdata3/title'>>>>> White Test'/rname={"Results:"}
                           /cname={"Chi-sq","df","p"}/format=!dcms.
                    print /title'Note   : Design = Intercept + IVs + square of IVs + interaction '+
                        'between IVs.'/space=0.
                else .
                    print /title'>>>>> White Test'/format=a8.
                    print /title'Note   : This test could be performed due to '+
                        'multicollinearity.'/space=0.
                end if.
            end if.
            compute tdata1=csum(tdata9).
            compute tdata3=tdata9/(tdata1/ncase).
            @candr crdata={tdata3,edtx2}.
            do if (candrn=0).
                @reg3 @rgy=tdata3/@rgx=edtx2/@rgci=eregci/@rg1vdr=0
                            /@rg1hc=9/@rgALL=1/@rg2cvg=lgconver.
                compute tdata1=csum(rg1e&**2).
                compute tdata4=make(ncase,ncase,1).
                compute tdata5=ident(ncase).
                compute tdata6=t(tdata3)*(tdata5-((1/ncase)*tdata4))*tdata3.
                compute bpregss=tdata6-tdata1.
                compute bpksq=1-(tdata1/tdata6).
                compute bptest=0.5*bpregss.
                compute bpsig=1-chicdf(bptest,nvarx2).
                compute bpktest=ncase*bpksq.
                compute bpksig=1-chicdf(bpktest,nvarx2).
                compute tdata1={'Chi-sq','df','p'}.
                print {bptest,nvarx2,bpsig}/title '>>>>> Breusch-Pagan Test'
                    /rname={"Results:"}/cname=tdata1/foramt=!dcms.
                print {bpktest,nvarx2,bpksig}/title '>>>>> Koenker Test'
                    /rname={"Results:"}/cname=tdata1/foramt=!dcms.
            else.
                print /title'>>>>> Breusch-Pagan Test'/format=a8.
                print /title'Note   : Unexpected error.'/space=0.
                print /title'>>>>> Koenker Test'/format=a8.
                print /title'Note   : Unexpected error.'/space=0.
            end if.
            compute tdata1=((csum(tdata8&**4)/ncase)/((csum(tdata8&**2)/ncase)**2)-3)*
                sqrt(ncase/24).
            compute tdata3=1-cdfnorm(abs(tdata1)).
            print {tdata1,tdata3}/title '>>>>> hhet Test (Klein et al., 2016)'/format=!dcms
                /rname={"Results:"}/cname={'Z','p'}.
            print /title'Note   : hhet Test is a one-tailed test.'/space=0.
        end if.
        compute tdata1=!mrplot1.
        compute tdata2=!mrplot2.
        compute tdata1=tdata1+tdata2.
        do if (tdata1>0 and dichy<>2).
            print /title=!plotitle.
            print /title=!plonote/space=0.
            do if (!mrplot1=1).
                print /title '====================> Histogram of Residual'/format=a8.
                    save tdata8/outfile='EReg_MR_res1.sav'/var=Residual.
                print /title 'get file="EReg_MR_res1.sav".'.
                print /title 'Frequencies variables=Residual/format=notable/histogram normal.'
                    /space=0.
            end if.
            do if (!mrplot2=1).
                print /title '====================> Predicted Value vs. Residual'/format=a8.
                    save {tdata7,tdata8}/outfile='EReg_MR_res2.sav'/var=Predicted Residual.
                print /title 'get file="EReg_MR_res2.sav".'.
                print /title 'Graph /scatterplot(bivar)=Predicted with Residual.'
                    /space=0.
            end if.
        end if.
        compute tdata1=!mrca+!mrda+!mrallsub.
        do if (tdata1>0 and nvarx2>1).
            compute mrnvar=nvarx2.
            compute mrnvar2=2**mrnvar-1.
            compute mrmodel=make(mrnvar2,(mrnvar+1),0).
            loop #i=1 to mrnvar2.
                compute tdata2=mrnvar2+1-#i.
                compute mrmodel(#i,1)=#i.
                loop #i2=1 to mrnvar.
                    compute tdata1=mod(tdata2,2).
                    compute tdata2=(tdata2-tdata1)/2.
                    compute mrmodel(#i,(#i2+1))=tdata1.
                end loop.
            end loop.
            compute allx2={make(ncase,1,1),edtx2}.
            compute nvarxy=nvarx2+1.
            do if (dichy<>2).
               compute allb=ginv(allx2)*edty.
               compute allymn=ncase*(csum(edty)/ncase)&**2.
               compute alltss=cssq(edty)-allymn.
               compute allrss=t(allb)*t(allx2)*edty-allymn.
               compute allmse=(alltss-allrss)/(ncase-nvarxy).
               compute allsubo=make(mrnvar2,9,0).
               loop #ab=1 to mrnvar2.
                  compute allxdt=make(ncase,1,1).
                  compute allnvar=1.
                  loop #ab2=1 to mrnvar.
                     do if (mrmodel(#ab,(#ab2+1))=1).
                        compute allxdt={allxdt,edtx2(:,#ab2)}.
                        compute allnvar=allnvar+1.
                     end if.
                  end loop.
                  compute allb2=ginv(allxdt)*edty.
                  compute allrss2=t(allb2)*t(allxdt)*edty-allymn.
                  compute alless2=alltss-allrss2.
                  compute allmse2=alless2/(ncase-allnvar).
                  compute allrsq2=allrss2/alltss.
                  compute allarsq2=1-(1-allrsq2)*(ncase-1)/(ncase-allnvar).
                  compute allcp=alless2/allmse+2*allnvar-ncase.
                  compute allaic=ncase*ln(alless2/ncase)+2*allnvar.
                  compute allbic=ncase*ln(alless2/ncase)+allnvar*ln(ncase).
                  compute allaicc=allaic+2*allnvar*(allnvar+1)/(ncase-allnvar-1).
                  compute allabic=ncase*ln(alless2/ncase)+ln((ncase+2)/24)*allnvar.
                  compute allsubo(#ab,1:9)={#ab,allrsq2,allarsq2,allmse2,
                        allaic,allbic,allcp,allaicc,allabic}.
               end loop.
               compute tdata1=!mrassot.
               compute tdata2=allsubo(:,(tdata1+1)).
               do if (tdata1<2.5).
                  @sortv sortv=tdata2/sorttp=-1.
                  else.
                  @sortv sortv=tdata2/sorttp=1.
               end if.
               compute mrmodel2=mrmodel.
               compute allsubo2=allsubo.
               loop #ab=1 to (nvarx2+1).
                  compute mrmodel2(sort2,#ab)=mrmodel(:,#ab).
               end loop.
               loop #ab=1 to 9.
                  compute allsubo2(sort2,#ab)=allsubo(:,#ab).
               end loop.
               compute tdata3={'R-square','Adj.R-sq','MSE','AIC','BIC','Cp'}.
               compute tdata2={'Sorting','by',tdata3(1,tdata1),'.'}.
               do if (!mrallsub=1).
                   print tdata2/title'********************************* All Subsets '+
                      'Regression ********************************'/rname={'Note   :'}/format=a8.
                   do if (nvarx2<8).
                      compute tdata1={'Model.No','R-square','Adj.R-sq','MSE','AIC','BIC','Cp',
    'AICc',
    'aBIC'}.
                      print allsubo2/title'====================> Fit'
                         /cname=tdata1/rname={'Results:'}/format=!dcms.
                      compute tdata1={'Model.No',nmx2}.
                 print mrmodel2/title'====================> Model and Variable (1 = Selected)'
                         /cname=tdata1/rname={'Results:'}/format=!dcms.
                      else.
                      save {mrmodel2,allsubo2(:,2:9)}/outfile='EReg_AllsubsetsR.sav' 
                            /var= Model !e_x Rsq adjRsq MSE AIC BIC Cp AICc aBIC.
                      print /title 'Note   : 1 = Selected.'/space=0.
                      print /title=!plonote.
                      print /title 'get file="EReg_AllsubsetsR.sav".'.
                   end if.
               end if.
               do if (!mrca=1).
                  compute cadata={mrmodel2(:,2:(nvarx2+1)),allsubo2(:,2)}.
                  compute catemp=rsum(cadata(:,1:nvarx2)).
                  compute catemp={catemp,mod(catemp,2)}.
                   loop #i=1 to mrnvar2.
                      compute catemp2=make(mrnvar2,1,0).
                      do if (catemp(#i,1)=nvarx2).
                         loop #i2 = 1 to mrnvar2.
                            do if (catemp(#i,2)=catemp(#i2,2)).
                               compute catemp2(#i2,1)=-1.
                               ELSE. 
                               compute catemp2(#i2,1)=1.
                            end if.
                         end loop.
                         else.
                         loop #i2=1 to mrnvar2.  
                               compute catemp3=1.
                               loop #i3=1 to nvarx2.
                                 do if (cadata(#i,#i3)=0).
                                    do if (cadata(#i2,#i3)=1).
                                       compute catemp3=catemp3*1.
                                       ELSE.
                                       compute catemp3=0.
                                    end if.
                                 end if.
                               end loop.
                               do if (catemp3=1).
                                  do if (catemp(#i,2)=catemp(#i2,2)).
                                     compute catemp2(#i2,1)=-1.
                                     ELSE. 
                                     compute catemp2(#i2,1)=1.
                                  end if.
                               end if .
                         end loop.
                      end if.
                      compute cadata(#i,(nvarx2+1))=csum(allsubo2(:,2)&*catemp2(:,1)).
                   end loop.
                   compute tdata1=mod(nvarx2,2).
                   do if (tdata1=1).
                      compute cadata(:,(nvarx2+1))=cadata(:,(nvarx2+1))*(-1).
                   end if.
                   @sortv sortv=catemp(:,1)/sorttp=1.
                   loop #i3=1 to (nvarx2+1).
                      compute cadata(sort2,#i3)=cadata(:,#i3).
                   end loop.
                   compute tdata1={nmx2(1,1:nvarx2),'R-square','%'}.
                   print /title='********************************** Commonality Analysis '+
                       '*********************************'/format=a8.
                   compute cacudata={cacudata,mrcr(2:(nvarx2+1),1)&**2-cacudata}.
                   print cacudata/cname={'Unique','Common'}/rname=nmx2/format=!dcms
                      /title'====================> Unique and Common R-square'.
                   compute cadata={cadata,(cadata(:,nvarx2+1))/mrrsq*100}.
                   do if (nvarx2<8).
                   print cadata/cname=tdata1/rname={'Results:'}/format=!dcms
                     /title'====================> Commonality Coefficient (1 = Selected)'.
                   else.
                   save cadata/outfile='EReg_Commonality.sav' /var= !e_x Rsquare Percent.
                   print /title'====================> Commonality Coefficient (1 = Selected)'.
                   print /title=!plonote/space=0.
                   print /title 'get file="EReg_Commonality.sav".'.
                   end if.
                end if.
                do if (!mrda=1).
                   compute dadata={mrmodel2(:,2:(nvarx2+1)),allsubo2(:,2)}.
                   compute dandata=ncol(dadata).
                   compute daprint1=make(nvarx2,(nvarx2+1),0).
                   compute daprint2=make(nvarx2,nvarx2,0).
                   compute danvar=rsum(dadata(:,1:nvarx2)).
                   loop #da= 1 to nvarx2.
                       compute datmp8=0.
                       do if (#da=1).
                           compute datmp1=make(2,dandata,0).
                           compute datmp3=2.
                           else.
                           compute datmp1=datmp2.
                           compute datmp3=nrow(datmp1).
                       end if.
                       compute datmp2=make(1,dandata,0).
                       compute datmp4=1.
                       loop #da2=1 to mrnvar2.
                           do if  (danvar(#da2,1)=#da).
                               compute datmp2={datmp2;dadata(#da2,:)}.
                               compute datmp4=datmp4+1.
                           end if.
                       end loop.
                       loop #da3=2 to datmp3.
                           compute datmp6={0,0}.
                           loop #da4=2 to datmp4.
                               compute datmp5=1.
                               loop #da5 =1 to nvarx2.
                                   do if (datmp1(#da3,#da5)=1).
                                       do if (datmp2(#da4,#da5)=1).
                                           compute datmp5=datmp5*1.
                                           else.
                                           compute datmp5=0.
                                       end if.
                                   end if.
                               end loop.
                               do if (datmp5=1).
                                   loop #da6=1 to nvarx2.
                                       do if (datmp1(#da3,#da6)=0 and datmp2(#da4,#da6)=1).
                                           compute datmp6={datmp6;#da6,
                                              (datmp2(#da4,(nvarx2+1))-datmp1(#da3,(nvarx2+1)))}.
                                       end if.
                                   end loop.
                               end if.
                           end loop.
                           compute datmp7=nrow(datmp6).
       loop #da7=2 to datmp7.
           loop #da8=2 to datmp7.
               do if (#da7<#da8).
                   do if (datmp6(#da7,2)>=datmp6(#da8,2)).
                       compute daprint2(datmp6(#da7,1),datmp6(#da8,1))=
                            {daprint2(datmp6(#da7,1),datmp6(#da8,1))+1}.
                       else.
                       compute daprint2(datmp6(#da8,1),datmp6(#da7,1))=
                            {daprint2(datmp6(#da8,1),datmp6(#da7,1))+1}.
                   end if.
               end if.
           end loop.
           compute daprint1(datmp6(#da7,1),(#da+1))=
                {daprint1(datmp6(#da7,1),(#da+1))+datmp6(#da7,2)}.
           compute datmp8=datmp8+1.
       end loop.
                       end loop.
                       compute daprint1(:,(#da+1))=(daprint1(:,(#da+1))/datmp8)*nvarx2.
                   end loop.
                   compute daprintd=daprint2(1,nvarx2)+daprint2(nvarx2,1).
                   compute daprint2=daprint2/daprintd.
                   compute daprint1(:,1)=rsum(daprint1(:,2:(nvarx2+1)))/(nvarx2).
                   print /title'*********************************** Dominance Analysis '+
                            '**********************************'.
compute tdata2={'G.D','K=0','K=1','K=2','K=3','K=4','K=5','K=6','K=7','K=8','K=9',
'K=10','K=11','K=12','K=13','K=14','K=15','K=16','K=17','K=18','K=19','K=20','K=21',
'K=22','K=23','K=24','K=25','K=26','K=27','K=28','K=29','K=30','K=31','K=32','K=33',
'K=34','K=35','K=36','K=37','K=38','K=39','K=40','K=41','K=42','K=43','K=44','K=45',
'K=46','K=47','K=48','K=49','K=50','K=51','K=52','K=53','K=54','K=55','K=56','K=57',
'K=58','K=59','K=60','K=61','K=62','K=63','K=64','K=65','K=66','K=67','K=68','K=69',
'K=70','K=71','K=72','K=73','K=74','K=75','K=76','K=77','K=78','K=79','K=80','K=81',
'K=82','K=83','K=84','K=85','K=86','K=87','K=88','K=89','K=90','K=91','K=92','K=93',
'K=94','K=95','K=96','K=97','K=98','K=99'}.
         print daprint1 /title '====================> General Dominance (G.D) and '+
                 'Conditional Dominance'   /rname=nmx2/cname=tdata2/format=!dcms.
         print /title'Note   : K is the number of IVs in model.'/space=0.
         compute tdata2=nmx2(1,1:nvarx2).
         print daprint2 /title '====================> Complete Dominance (row vs. column)'
            /cname=tdata2/rname=nmx2/format=!dcms.
                end if.
            end if.
            do if (dichy=2).
                compute allsubo=make(mrnvar2,7,0).
                loop #ab=1 to mrnvar2.
                    compute allxdt=make(ncase,1,1).
                    compute allxnvar=1.
                    loop #ab2=1 to mrnvar.
                        do if (mrmodel(#ab,(#ab2+1))=1).
                            compute allxdt={allxdt,edtx2(:,#ab2)}.
                            compute allxnvar=allxnvar+1.
                        end if.
                    end loop.
                    compute allxdt=allxdt(:,2:allxnvar).
                    @reg3 @rgy=edty/@rgx=allxdt/@rgci=eregci/@rg1vdr=0
                                /@rg1hc=9/@rgALL=0/@rg2cvg=lgconver.
                    compute allchisq=rg2hstr(1,2)-rg2LL2.
                    compute allaic= rg2LL2+2*(rgnvar+1).
                    compute allbic= rg2LL2+rgnvar*ln(ncase).
                    compute allaicc= allaic+2*(rgnvar+1)*(rgnvar+2)/(ncase-rgnvar-2).
                    compute allabic= rg2LL2+rgnvar*ln((ncase+2)/24).
                    compute allsubo(#ab,1:7)={#ab,rg2LL2,allchisq,allaic,allbic,allaicc,allabic}.
                end loop.
                compute tdata1=!mrassot.
                do if (tdata1>5).
                    compute tdata1=5.
                end if.
                do if (tdata1<4).
                    compute tdata1=2.
                end if.
                compute tdata2=allsubo(:,tdata1).
                @sortv sortv=tdata2/sorttp=1.
                compute mrmodel2=mrmodel.
                compute allsubo2=allsubo.
                loop #ab=1 to (nvarx2+1).
                     compute mrmodel2(sort2,#ab)=mrmodel(:,#ab).
                end loop.
                loop #ab=1 to 7.
                     compute allsubo2(sort2,#ab)=allsubo(:,#ab).
                end loop.
                compute tdata3={'-2LL','-2LL','AIC','BIC','AICc','aBIC'}.
                compute tdata2={'Sorting','by',tdata3(1,(tdata1-1)),'.'}.
                do if (!mrallsub=1).
                    print tdata2/title'********************************* All Subsets '+
                         'Regression ********************************'/rname={'Note   '+
    ':'}/format=a8.    
     do if (nvarx2<8).
          compute tdata1={'Model.No','-2LL','Chi-sq','AIC','BIC','AICc','aBIC'}.
          print allsubo2/title'====================> Fit'
             /cname=tdata1/rname={'Results:'}/format=!dcms.
          compute tdata1={'Model.No',nmx2}.
          print mrmodel2/title'====================> Model and Variable (1 = Selected)'
             /cname=tdata1/rname={'Results:'}/format=!dcms.
          else.
          save {mrmodel2,allsubo2(:,2:7)}/outfile='EReg_AllsubsetsR.sav' 
                /var= Model !e_x Negtive2LL Chisq AIC BIC AICc aBIC.
          print /title 'Note   : 1 = Selected.'/space=0.
          print /title=!plonote.
          print /title 'get file="EReg_AllsubsetsR.sav".'.
     end if.
                end if.
            end if.
        end if.
    end if.
end if.
do if (emod=3).
    !let !tt=!concat ('Focal IV (X): ',!head(!e_x))
    !let !tt = !quote(!tt).
    print /title=!tt.
    !let !tt=!concat('Quadratic term : X^2',' ')
    !let !tt = !quote(!tt).
    print /title=!tt/format=a8/space=0.
    compute qedttp=!qedttp.
    compute qepap=!qepap.
    compute simslp=qepap.
    compute qejn=!qejn.
    compute qejnp=!qejnp.
    compute qesave=!qesave.
    compute qerbst=!qerbst.
    do if (dichy<>2).
        do if (qerbst<>9).
            !let !tt = !concat('Standard errors for regression coefficients : HC',!qerbst)
            !let !tt = !quote (!tt)
            print /title=!tt/format=a8.
        end if .
    end if.
    do if (qedttp=1).
        print /title'Data type : Raw'.
    else if (qedttp=2).
        !let !tt =!concat('Data type : Mean center (only for ',!head(!e_x),')')
        !let !tt = !quote(!tt).
        print /title=!tt.
    else if (qedttp=3).
        print /title 'Data type : Standardization (all variables) '.
    end if.
    do if (nvarx2>(ncase-3)).
        compute ero={ero;7}.
    end if.
    compute sero=csum(ero).
    do if (sero=0).
        compute tdata1={edty,edtx2}.
        compute qevcv=ncase2*(sscp(tdata1)-((t(csum(tdata1))*csum(tdata1))/ncase)).
        compute qesd=diag(qevcv)&**0.5.
        compute qemean=csum(tdata1)/ncase.
        do if (qedttp=2).
            compute edtx1=edtx1-csum(edtx1)/ncase.
        else if (qedttp=3).
            loop #i=1 to (nvarx2+1).
                compute tdata2=tdata1(:,#i).
                compute tdata1(:,#i)=(tdata1(:,#i)-qemean(1,#i))/qesd(#i,1).
            end loop.
        end if.
        compute edty=tdata1(:,1).
        compute edtx1=tdata1(:,2).
        compute edtx2=tdata1(:,2:(nvarx2+1)).
        do if (nvarx2>1).
            compute edtco=edtx2(:,2:nvarx2).
        end if.
        compute edtx12=edtx1&**2.
        compute qeregx={edtx1,edtx12}.
        do if (nvarx2>1).
            compute qeregx={qeregx,edtco}.
        end if.
        compute tdata1={edty,qeregx}.
        @candr crdata=tdata1.
        compute ero={ero;candrn}.
        compute sero=csum(ero).
        do if (sero=0).
            compute tdata1={edty,edtx2}.
            compute qevcv=ncase2*(sscp(tdata1)-((t(csum(tdata1))*csum(tdata1))/ncase)).
            compute qesd=sqrt(diag(qevcv)).
            compute qemean=csum(tdata1)/ncase.
            compute qemax=cmax(tdata1).
            compute qemin=cmin(tdata1).
            compute tdatanm={'Constant',nmx1,'X^2'}.
            do if (nvarx2>1).
                compute tdatanm={tdatanm,nmx2(1,2:(nvarx2))}.
            end if.
            compute qedes={t(qemean),qesd,t(qemin),t(qemax)}.
            compute tdatanm2={nmy,nmx2,'       .'}.
            print qedes/title='*************************'+
                '******** Descriptive Statistics ********************************'
                /rname=tdatanm2/clabels='Mean' 'Std.Dev'  'Min' 'Max' /format=!dcms.
            print /title=!rgline.
            print /title=!rgytitle.
            @reg3 @rgy=edty/@rgx=qeregx/@rgci=eregci/@rg1vdr=1/@rg1hc=qerbst
                /@rgALL=1/@rg2cvg=lgconver.
            do if (dichy=2).
                print rgycode/title=!rgycode/format=!dcms/rname=prtrslt/cname=rgycodep.
                print rg2ifm/title '>>>>> Model Summary'/format=!dcms/rname=prtrslt/cname=rg2ifmp.
                print rg2class/title '>>>>> Classification'/format=!dcms
                    /rname=rg2clspr/cname=rg2clspc.
                compute qtout=ppnd16.
                @polystdb @polyb=rg2coef(:,8)/@polysd=rg2sd(4:(4+rgnvar),:)/@polym=11.
                @polystdb @polyb=rg2coef(:,9)/@polysd=rg2sd(4:(4+rgnvar),:)/@polym=11.
                print rg2coef/title '>>>>> Coefficients'/format=!dcms
                    /rname=tdatanm/cname=rgbp.
            else.
                print rg1ifm/title '>>>>> Model Summary'/format=!dcms/rname=prtrslt/cname=rg1ifmp.
                @polystdb @polyb=rg1coef(:,7)/@polysd=rg1sd/@polym=11.
                print rg1coef/title '>>>>> Coefficients'/format=!dcms
                    /rname=tdatanm/cname=rgbp.
            end if.
            compute fg1=rgb(3,1).
            compute fb1=rgb(2,1).
            /** Utest **/.
            compute fs11=rgbcov(2,2).
            compute fs12=rgbcov(2,3).
            compute fs22=rgbcov(3,3).
            print /title'************************ Test of (Inverted) U-shaped Relationship '+
                '***********************'.
            compute fd=(fs12**2-fs22*fs11)*(qtout**2)+fg1**2*fs11+fb1**2*fs22-2*fs12*fb1*fg1.
            compute extp={-rgb(2,1)/rgb(3,1)*0.5,-999,-999}.
            do if (fd>0).
                compute extp(1,2)=0.5*(fs12*qtout**2-fb1*fg1-qtout*sqrt(fd))
                    /(fg1**2-fs22*qtout**2).
                compute extp(1,3)=0.5*(fs12*qtout**2-fb1*fg1+qtout*sqrt(fd))
                    /(fg1**2-fs22*qtout**2).
            end if.
            compute tdata1=sqrt((fs11*fg1**2-2*fb1*fg1*fs12+fs22*fb1**2)/(4*fg1**4)).
            compute extp={extp;extp(1,1),(extp(1,1)-tdata1*qtout),(extp(1,1)+tdata1*qtout)}.
            print extp/title'====================> Turning Point'/rname={'Fieller:','Delta  :'}
                    /cname={'Point.X','LLCI','ULCI'}/format=!dcms.
            /** Simple slope **/.
            compute qesbcov={rgbcov(2,2),rgbcov(2,3);rgbcov(3,2),rgbcov(3,3)}.
            compute qexmax=qedes(2,4).
            compute qexmin=qedes(2,3).
            compute qexincn=21.
            compute qexinc=(qexmax-qexmin)/(qexincn-1).
            loop #jn=1 to qexincn.
                do if (#jn=1).
                    compute qexjnlp=qexmin.
                    else.
                    compute qexjnlp={qexjnlp;(qexmin+qexinc*(#jn-1))}.
                end if.
            end loop.
            /** jn critical value **/.
            compute fd2=(4*qtout**2*fs12-4*fb1*fg1)**2-4*(4*qtout**2*fs22-4*fg1**2)*
                                (qtout**2*fs11-fb1**2).
            do if (fd2=0). 
                compute qejncd=(4*fb1*fg1-4*qtout**2*fs12)/(2*(4*qtout**2*fs22-4*fg1**2)).
            else if (fd2>0).
                compute tdata1=(4*fb1*fg1-4*qtout**2*fs12-fd2**0.5)
                    /(2*(4*qtout**2*fs22-4*fg1**2)).
                compute tdata2=(4*fb1*fg1-4*qtout**2*fs12+fd2**0.5)
                    /(2*(4*qtout**2*fs22-4*fg1**2)).
                compute qejncd={tdata1;tdata2}.
            end if.
            do if (fd2>=0).
                compute qexjnlp={qexjnlp;qejncd}.
                compute tdata1=nrow(qexjnlp).
                @sortv sortv=qexjnlp/sorttp=1.
                compute tdata2=0.
                loop #i= 1 to tdata1.
                    do if (qexjnlp(#i,1)>=qexmin and qexjnlp(#i,1)<=qexmax).
                        compute tdata2={tdata2;qexjnlp(#i,1)}.
                    end if.
                end loop.
                compute qejnn=nrow(tdata2)-1.
                compute qexjnlp=tdata2(2:(qejnn+1),1).
            else.
                compute qejnn=qexincn.
            end if.
            do if (qepap=3).
                compute qepapv={!qepapv}.
                compute qeslpp={t(qepapv)}.
                compute qepaprnm={'Custom.V'}.
            else if (qepap=1).
                compute qeslpp={(qedes(2,1)-qedes(2,2));qedes(2,1);(qedes(2,1)+qedes(2,2))}.
            else if (qepap=2).
                compute tdata1=edtx1.
                @sortv sortv=tdata1/sorttp=1.
                @prcn prcny=tdata1/prcnp=16.
                compute qeslpp=percnv.
                @prcn prcny=tdata1/prcnp=50.
                compute qeslpp={qeslpp;percnv}.
                @prcn prcny=tdata1/prcnp=84.
                compute qeslpp={qeslpp;percnv}.
            end if.
            compute qexjnlp={qexjnlp;qeslpp}.
            compute qejnn2=nrow(qexjnlp).
            compute fxjn=make(qejnn2,7,-999).
            compute fxjn(:,1)=qexjnlp.
            compute fxjn(:,2)=fxjn(:,1)*2*rgb(3,1)+rgb(2,1).
            compute qesbdata={make(qejnn2,1,1),fxjn(:,1)*2}.
            compute fxjn(:,3)=sqrt(diag(qesbdata*qesbcov*t(qesbdata))).
            compute aa=fxjn(:,2)&/fxjn(:,3).
            compute fxjn(:,4)=fxjn(:,2)&/fxjn(:,3).
            do if (dichy=2).
                compute fxjn(:,5)=2*(1-cdfnorm(abs(fxjn(:,4)))).
            else.
                compute fxjn(:,5)=2*(1-tcdf(abs(fxjn(:,4)),rg1df2)).
            end if.
            compute fxjn(:,6:7)={(fxjn(:,2)-(fxjn(:,3)*qtout)),(fxjn(:,2)+(fxjn(:,3)*qtout))}.
            compute slpend={fxjn(1,:);fxjn(qejnn,:)}.
            compute slpend={slpend,slpend(:,5)/2}.
            compute slpnm={nmx1,'Slope','SE','t','p','LLCI','ULCI','Overall'}.
            do if (dichy=2).
                compute slpnm(1,4)={'Z'}.
            end if.
            print slpend/title'====================> Slopes'
                /rname={'Min.End:','Max.End:'}/cname=slpnm/foramt=!dcms.
            print /title'************************************* Simple Slope '+
                '**************************************'.
            compute sslpcnm=slpnm(1,1:7).
            compute sslprnm={make(23,1,' ');'       .'}.
            compute qepapd=fxjn((qejnn+1):qejnn2,:).
            compute fxjn=fxjn(1:qejnn,:).
            print qepapd/title '====================> Pick A Point '
                /cname=sslpcnm/rname=sslprnm/foramt=!dcms.
            do if (qejn=1).
                 PRINT fxjn/title '====================> Johnson-Neyman'
                         /cname=sslpcnm/rname=sslprnm/format=!dcms.
             end if.
             do if (qejnp=1).
                PRINT /title '====================> Code for Plot (Johnson-Neyman)'.
                print /title=!plonote/space=0.
                !let !tt=!e_x
                !let !tt2=!head(!tt)
                !let !tt = !concat('Data list free/',!tt2,!BLANKS(1),'Slope LLCI ULCI .')
                !let !tt = !quote(!tt)
                print /title=!tt /foramt=a8.
                print {fxjn(:,1:2),fxjn(:,6:7)}/title 'Begin data.'/format=!dcms/space=0.
                print /title 'End data.'/format=a8/space=0.
                !let !tt = !concat('Graph/scatterplot(overlay)=',
                !tt2,!blanks(1),!tt2,!blanks(1),!tt2,!blanks(1),' With ','Slope LLCI ULCI (Pair).')
                !let !tt = !quote(!tt)
                print /title=!tt/format=a8/space=0.
             end if.
             do if (qesave=1).
                print /title '*************************************** Save Data '+
                    '***************************************'.
                save {eregid,edty,edtx2,edtx2(:,1)&**2}/outfile=*/vars=EReg_ID !e_y !e_x.
                compute prnm={colnm(1,(nvarx2+3)),'X^2'}.
                print prnm/title'>>>>> Column Names vs. Variables' /cname={'COL '+
                    'Name','Variable'}  /rname={'List   :'}/format=a8.
             end if.
        end if.
    end if.
end if.
do if (emod=4).
    do if (dichy<>2).
        compute modrbst=!modrbst.
        do if (modrbst<>9).
            !let !tt = !concat('Standard errors for regression coefficients : HC',!modrbst)
            !let !tt = !quote (!tt)
            print /title=!tt/format=a8.
        end if .
    end if.
    compute tdata1=nmx2(1,1).
    print tdata1/title 'Focal IV (X) is '/format=a8.
    compute modm=!modm.
    do if (nvarx2=1).
        print /title 'Moderator does not exist. '/format=a8/space=0.
        compute ero={ero;14}.
    else.
        compute tdata2={nmx2(1,2)}.
        print tdata2/title 'Moderator (M) is '/format=a8/space=0.
        compute intnm={'X^2','=',tdata1,'*',tdata1,' ',' ';'M^2','=',tdata2,'*',tdata2,' ',' ';
        'X*M','=',tdata1,'*',tdata2,' ',' ';'X^2*M','=',tdata1,'*',tdata1,'*',
        tdata2;'X*M^2','=',tdata1,'*',tdata2,'*',tdata2}.
        do if (modm=6).
            compute intnmp=intnm(3,:).
            else if (modm=7).
            compute intnmp={intnm(3,:);intnm(1,:)}.
            else if (modm=8).
            compute intnmp={intnm(3,:);intnm(1,:);intnm(4,:)}.
            else if (modm=9).
            compute intnmp={intnm(3,:);intnm(2,:)}.
            else if (modm>=10).
            compute intnmp={intnm(3,:);intnm(2,:);intnm(5,:)}.
        end if.
        print intnmp/title 'Product term(s) key:'/rname={'List   :'}/format=a8/space=0.
    end if.
    do if (!moddttp=1).
        print /title'Data type : Raw'.
        else if (!moddttp=2).
        print /title'Data type : Mean center (only for X and M)'.
        else if (!moddttp=3).
        print /title 'Data type : Standardization (all variables) '.
    end if.
    do if (modm=6).
        compute tdata1=nvarx2+1.
    else if (modm=7).
        compute tdata1=nvarx2+2.
    else if (modm=8).
        compute tdata1=nvarx2+3.
    else if (modm=9).
        compute tdata1=nvarx2+2.
    else if (modm=10).
        compute tdata1=nvarx2+3.
    end if.
    do if (tdata1>(ncase-2)).
        compute ero={ero;7}.
    end if.
    compute sero=csum(ero).
    do if (sero=0).
        compute modx=edtx2(:,1).
        compute modmo=edtx2(:,2).
        compute tdata1={edty,edtx2}.
        compute modvcv=ncase2*(sscp(tdata1)-((t(csum(tdata1))*csum(tdata1))/ncase)).
        compute modsd=sqrt(diag(modvcv)).
        compute modd=inv(mdiag(modsd)).
        compute modcr=modd*modvcv*modd.
        compute modmean=t(csum(tdata1)/ncase).
        compute modpap=!modpap.
        compute simslp=modpap.
        do if (nvarx2>2).
            compute modcov=edtx2(:,3:nvarx2).
        end if.
        do if (!moddttp=2).
            compute modx=modx-make(ncase,1,modmean(2,1)).
            compute modmo=modmo-make(ncase,1,modmean(3,1)).
            compute tdata1(:,2:3)={modx,modmo}.
            else if (!moddttp=3).
            loop #mod=1 to (nvarx2+1).
                compute tdata1(:,#mod)=(tdata1(:,#mod)-
                    make(ncase,1,modmean(#mod,1)))/modsd(#mod,1).
            end loop.
            compute edty=tdata1(:,1).
            compute modx=tdata1(:,2).
            compute modmo=tdata1(:,3).
            do if (nvarx2>2).
                compute modcov=tdata1(:,4:(nvarx2+1)).
            end if.
        end if.
        compute modmean=t(csum(tdata1)/ncase).
        do if (!moddttp=3).
            compute modsd=make((nvarx2+1),1,1).
        end if.
        compute modmin=t(cmin(tdata1)).
        compute modmax=t(cmax(tdata1)).
        compute moddes={modmean,modsd,modmin,modmax}.
        compute tdatanm2={nmy,nmx2,'       .'}.
        print moddes/title='*************************'+
            '******** Descriptive Statistics ********************************'
            /rname=tdatanm2/clabels='Mean' 'Std.Dev'  'Min' 'Max' /format=!dcms.
        release tdata1.
        compute modinta={modx&**2,modmo,modmo&**2,
            modx&*modmo,modx&**2&*modmo,modx&*modmo&**2}.
        compute modintn={'X^2',nmx2(1,2),'M^2','X*M','X^2*M','X*M^2'}.
        compute modints={0,1,0,1,0,0;1,1,0,1,0,0;1,1,0,1,1,0;0,1,1,1,0,0;0,1,1,1,0,1}.
        compute modints=modints((modm-5),:).
        compute modxnm={'Constant',nmx1}.
        loop #mod=1 to 6.
            do if (modints(1,#mod)=1).
                compute modx={modx,modinta(:,#mod)}.
                compute modxnm={modxnm,modintn(1,#mod)}.
            end if.
        end loop.
        do if (nvarx2>2).
            compute modx={modx,modcov}.
            compute modxnm={modxnm,nmx2(1,3:nvarx2)}.
        end if.
        print /title=!rgline.
        print /title=!rgytitle/format=a8.
        @reg3 @rgy=edty/@rgx=modx/@rgci=eregci/@rg1vdr=1/@rg1hc=!modrbst
                /@rgALL=1/@rg2cvg=lgconver.
        do if (dichy=2).
            compute qtout=ppnd16.
            print rgycode/title=!rgycode/format=!dcms/rname=prtrslt/cname=rgycodep.
            print rg2ifm/title '>>>>> Model Summary'/format=!dcms/rname=prtrslt/cname=rg2ifmp.
            print rg2class/title '>>>>> '+
                    'Classification'/format=!dcms/rname=rg2clspr/cname=rg2clspc.
            @polystdb @polyb=rg2coef(:,8)/@polysd=rg2sd(4:(4+rgnvar),:)/@polym=modm.
            @polystdb @polyb=rg2coef(:,9)/@polysd=rg2sd(4:(4+rgnvar),:)/@polym=modm.
            print rg2coef/title '>>>>> Coefficients'/format=!dcms
                /rname=modxnm/cname=rgbp.
            else.
            print rg1ifm/title '>>>>> Model Summary'/format=!dcms/rname=prtrslt/cname=rg1ifmp.
            @polystdb @polyb=rg1coef(:,7)/@polysd=rg1sd/@polym=modm.
            print rg1coef/title '>>>>> Coefficients'/format=!dcms
                /rname=modxnm/cname=rgbp.
        end if.
        do if (modpap=1).
            compute modpapv={moddes(3,1)-moddes(3,2);moddes(3,1);moddes(3,1)+moddes(3,2)}.
            compute modpapr=1.
          else if (modpap=2).
            @sortv sortv=modmo/sorttp=1.
            @prcn prcny=modmo/prcnp=16.
            compute modpapv=percnv.
            @prcn prcny=modmo/prcnp=50.
            compute modpapv={modpapv;percnv}.
            @prcn prcny=modmo/prcnp=84.
            compute modpapv={modpapv;percnv}.
            compute modpapr=1.
          else if (modpap=3).
            compute modpapv={!modpapv}.
            compute modpapr=nrow(modpapv).
            do if (modm=7 or modm=8).
                do if (modpapr>1).
                    compute modpapv=modpapv(1:2,1:2).
                    else.
                    compute modpapv=t(modpapv).
                    @sortv sortv=modpapv/sorttp=1.
                end if.
                else.
                compute modpapv=t(modpapv(1,:)).
                @sortv sortv=modpapv/sorttp=1.
                compute modpapr=1.
            end if.
        end if.
        compute modjnn=!modjnn.
        compute modmomin=modmin(3,1).
        compute modmomax=modmax(3,1).
        compute modmoinc=(modmomax-modmomin)/modjnn.
        compute modssn=modmomin.
        compute modxmin=modmin(2,1).
        compute modxmax=modmax(2,1).
        compute modxinc=(modxmax-modxmin)/modjnn.
        compute modqssn=modxmin.
        loop #ss=1 to modjnn.
            compute modssn={modssn;(modmomin+#ss*modmoinc)}.
            compute modqssn={modqssn;(modxmin+#ss*modxinc)}.
        end loop.
        compute modssnx=modssn.
        print /title'************************************* Simple Slope '+
            '**************************************'.
        compute modssnm={nmx2(1,2),'Slope','SE','t','p','LLCI','ULCI'}.
        do if (dichy=2).
            compute modssnm(1,4)={'Z'}.
        end if.
        compute modssrnm=make((modjnn+nrow(modpapv)+99),1,'        ').
        compute modssrnm={modssrnm;'       .'}.
        compute modssnm3={'Slope A','Slope B','A - B',modssnm(1,3:7)}.
        compute qtout2=qtout**2.
        /** simple slope **/.
        do if (modm=6).
            compute mods11=rgbcov(2,2).
            compute mods13=rgbcov(4,2).
            compute mods33=rgbcov(4,4).
            compute modb1=rgb(2,1).
            compute modb3=rgb(4,1).
            @eq.slv @eq.num={(qtout2*mods33-modb3**2),(2*qtout2*mods13-2*modb1*modb3),
                    (qtout2*mods11-modb1**2)}.
            do if (@eq.rtc=1).
                compute modssn={modssn;@eq.rt}.
                @num.out @no=modssn/@nmm={modmomin;modmomax}.
                @sortv sortv=modssn/sorttp=1.
            end if.
            compute modjnn=nrow(modssn).
            compute modssn={modssn;modpapv}.
            compute modjnn2=nrow(modssn).
            compute modssn={modssn,make(modjnn2,6,0)}.
            compute modb2={rgb(2,1);rgb(4,1)}.
            compute modbcov2={rgbcov(2,2),rgbcov(4,2);rgbcov(2,4),rgbcov(4,4)}.
            compute tdata1={make(modjnn2,1,1),modssn(:,1)}.
            compute modssn(:,2)=tdata1*modb2.
            compute modssn(:,3)=sqrt(diag(tdata1*modbcov2*t(tdata1))).
            compute modssn(:,4)=modssn(:,2)&/modssn(:,3).
            do if (dichy=2).
                compute modssn(:,5)=2*(1-cdfnorm(abs(modssn(:,4)))).
                else.
                compute modssn(:,5)=2*(1-tcdf(abs(modssn(:,4)),rg1df2)).
            end if.
            compute modssn(:,6:7)={(modssn(:,2)-modssn(:,3)*qtout),
                (modssn(:,2)+modssn(:,3)*qtout)}.
            compute modssn2=modssn((modjnn+1):modjnn2,:).
            compute modssn=modssn(1:modjnn,:).
            print modssn2/title '====================> Pick A Point '
                /cname=modssnm/rname=modssrnm/foramt=!dcms.
            do if (!modjn=1).
                print modssn/title '====================> Johnson-Neyman'
                        /cname=modssnm/rname=modssrnm/format=!dcms.
            end if.
  do if (!modjnp=1).
      print /title '====================> Code for Plot (Johnson-Neyman)'.
      compute modjnptt={'Data ','list ','free/',nmx2(1,2),' Slope ','LLCI','ULCI.'}.
      print /title=!plonote/format=a8/space=0.
      print modjnptt/title='  '/format=a8/space=0.
      print {modssn(:,1:2),modssn(:,6:7)}/title 'Begin data.'/format=!dcms/space=0.
      print /title 'End data.'/format=a8/space=0.
      compute tdata1={nmx2(1,2),nmx2(1,2),nmx2(1,2),'With';'Slope',' LLCI','ULCI','(Pair).'}.
      print tdata1/title="Graph/scatterplot(overlay)="/format=a8/space=0.
  end if.
        else if (modm=7 or modm=8).
            compute modpapvn=nrow(modpapv).
            compute modb=rgb(2:5,1)}.
            compute modbcov=rgbcov(2:5,2:5).
            do if (modm=8).
                compute modb=rgb(2:6,1).
                compute modbcov=rgbcov(2:6,2:6).
            end if.
            do if (modpapr=1).
                compute modpapv={make(modpapvn,1,modmean(2,1)),modpapv}.
            end if.
            compute modxssv={make(modpapvn,1,1),modpapv(:,1)*2,
                make(modpapvn,1,0),modpapv(:,2)}.
            do if (modm=8).
                compute modxssv={modxssv,modpapv(:,1)*2&*modpapv(:,2)}.
            end if.
            loop #md=1 to modpapvn.
                loop #md2=1 to modpapvn.
                    do if (#md<#md2).
                         compute modxssv={modxssv;modxssv(#md,:)-modxssv(#md2,:)}.
                    end if.
                end loop.
            end loop.
            compute modpapvm=modpapvn*(modpapvn-1)/2.
            compute modssn=make((modpapvn+modpapvm),8,-99).
            compute modssn(1:modpapvn,1:2)=modpapv.
            compute modssn(:,3)=modxssv*modb.
            compute tdata1=modpapvn.
            loop #md=1 to modpapvn.
                loop #md2=1 to modpapvn.
                    do if (#md<#md2).
                        compute tdata1=tdata1+1.
                        compute modssn(tdata1,1:2)={modssn(#md,3),modssn(#md2,3)}.
                    end if.
                end loop.
            end loop.
            compute modssn(:,4)=sqrt(diag(modxssv*modbcov*t(modxssv))).
            compute modssn(:,5)=modssn(:,3)&/modssn(:,4).
            do if (dichy=2).
                compute modssn(:,6)=2*(1-cdfnorm(abs(modssn(:,5)))).
                else.
                compute modssn(:,6)=2*(1-tcdf(abs(modssn(:,5)),rg1df2)).
            end if.
            compute modssn(:,7:8)={(modssn(:,3)-modssn(:,4)*qtout),
                (modssn(:,3)+modssn(:,4)*qtout)}.
            compute modssnm2={nmx2(1,1),modssnm}.
            print modssn(1:modpapvn,:)/title '====================> Pick A Point '
                    /cname=modssnm2/rname=modssrnm/foramt=!dcms.
            do if (modpapr=1).
                print /title 'Note   : The focal IV is evaluated at mean value.' /space=0.
                else.
                print modssn((modpapvn+1):(modpapvn+modpapvm),:)
                    /title '====================> Slope Difference Test'
                    /cname=modssnm3/rname=modssrnm/format=!dcms.
            end if.
            do if (modpapr=1).
                compute modpapv=modpapv(:,2).
                compute modssnm={nmx2(1,2),nmx2(1,1),modssnm(1,2:7)}.
                compute modjnptt={'Data ','list ','free/',nmx2(1,2),nmx2(1,1),' Slope '+
    '','LLCI','ULCI.'}.
                compute modjnpt2={nmx2(1,1),nmx2(1,1),nmx2(1,1),'With';'Slope',' '+
    'LLCI','ULCI','(Pair).'}.
                loop #md=1 to modpapvn.
                    compute modtm=modpapv(#md,1).
                    compute modjna=4*qtout2*modbcov(2,2)-4*modb(2,1)**2.
                    compute modjnb=4*(qtout2*modbcov(1,2)+qtout2*modtm*modbcov(2,4)-
                            modb(1,1)*modb(2,1)-modb(2,1)*modb(4,1)*modtm).
                    compute modjnc=qtout2*modbcov(1,1)+2*qtout2*modtm*modbcov(1,4)+
                            qtout2*modtm**2*modbcov(4,4)-modb(1,1)**2-2*modb(1,1)*modb(4,1)*modtm-
                            modb(4,1)**2*modtm**2.
                    do if (modm=8).
                        compute modjna=modjna+8*qtout2*modtm*modbcov(2,5)+
                                4*qtout2*modtm**2*modbcov(5,5)-
                                8*modb(2,1)*modb(5,1)*modtm-4*modb(5,1)**2*modtm**2.
                        compute modjnb=modjnb+4*(qtout2*modtm*modbcov(1,5)+
                                qtout2*modtm**2*modbcov(4,5)-modb(1,1)*modb(5,1)*modtm-
                                modb(4,1)*modb(5,1)*modtm**2).
                    end if.
                    compute modqssn2=modqssn.
                    @eq.slv @eq.num={modjna,modjnb,modjnc}.
                    do if (@eq.rtc=1).
                        compute modqssn2={modqssn2;@eq.rt}.
                        @num.out @no=modqssn2/@nmm={modxmin;modxmax}.
                        @sortv sortv=modqssn2/sorttp=1.
                    end if.
                    compute modjnn=nrow(modqssn2).
                    compute modqssn2={make(modjnn,1,1),modqssn2*2,
                        make(modjnn,1,0),make(modjnn,1,modtm)}.
                    do if (modm=8).
                        compute modqssn2={modqssn2,modqssn2(:,2)*modtm}.
                    end if.
                    compute modqssnp=modqssn2*modb.
                    compute modqssnp={modqssn2(:,2)/2,modqssnp,
                        sqrt(diag(modqssn2*modbcov*t(modqssn2)))}.
                    compute modqssnp={modqssnp,modqssnp(:,2)/modqssnp(:,3)}.
                    do if (dichy=2).
                        compute modqssnp={modqssnp,2*(1-cdfnorm(abs(modqssnp(:,4))))}.
                        else.
                        compute modqssnp={modqssnp,2*(1-tcdf(abs(modqssnp(:,4)),rg1df2))}.
                    end if.
                    compute modqssnp={modqssnp,(modqssnp(:,2)-modqssnp(:,3)*qtout),
                                                    (modqssnp(:,2)+modqssnp(:,3)*qtout)}.
                    compute modqssnp={make(modjnn,1,modtm),modqssnp}.
                    do if (!modjn=1).
                        print modqssnp/title '====================> Johnson-Neyman'
                                /cname=modssnm/rname=modssrnm/format=!dcms.
                    end if.
     do if (!modjnp=1).
         print /title '====================> Code for Plot (Johnson-Neyman)'.
         print /title=!plonote/format=a8/space=0.
         print modjnptt/title='  '/format=a8/space=0.
         print {modqssnp(:,1:3),modqssnp(:,7:8)}/title 'Begin data.'/format=!dcms/space=0.
         print /title 'End data.'/format=a8/space=0.
         print modjnpt2/title="Graph/scatterplot(overlay)="/format=a8/space=0.
     end if.
                end loop.
            end if.
            do if (!modjn2=1).
               compute modjnx={!modjnx}.
               compute modjnx=modjnx(1,:).
               compute modjnxn=ncol(modjnx).
               compute modssn2=modssnx.
               compute modssnm(1,1:2)=nmx2(1,1:2).
               compute modjnptt={'Data ','list ','free/',nmx2(1,1),
                    nmx2(1,2),' Slope ','LLCI','ULCI.'}.
               compute modjnpt2={nmx2(1,2),nmx2(1,2),
                    nmx2(1,2),'With';'Slope',' LLCI','ULCI','(Pair).'}.
               loop #jnx=1 to modjnxn.
                   compute modtm=modjnx(1,#jnx).
                   compute modjna=qtout2*modbcov(4,4)-modb(4,1)**2.
                   compute modjnb=2*qtout2*modbcov(1,4)+4*qtout2*modtm*modbcov(2,4)-
                                       2*modb(1,1)*modb(4,1)-4*modb(2,1)*modb(4,1)*modtm.
                   compute modjnc=qtout2*modbcov(1,1)+4*qtout2*modtm*modbcov(1,2)+
                                       4*qtout2*modtm**2*modbcov(2,2)-modb(1,1)**2-
                                       4*modb(1,1)*modb(2,1)*modtm-4*modb(2,1)**2*modtm**2.
                   do if (modm=8).
                       compute modjna=modjna+4*qtout2*modtm**2*modbcov(5,5)+
                                   4*qtout2*modtm*modbcov(4,5)-4*modb(4,1)*modb(5,1)*modtm-
                                   4*modb(5,1)**2*modtm**2.
                       compute modjnb=modjnb+4*qtout2*modtm*modbcov(1,5)+
                                   8*qtout2*modtm**2*modbcov(2,5)-4*modb(1,1)*modb(5,1)*modtm-
                                   8*modb(2,1)*modb(5,1)*modtm**2.
                   end if.
                   @eq.slv @eq.num={modjna,modjnb,modjnc}.
                   do if (@eq.rtc=1).
                       compute modssn2={modssn2;@eq.rt}.
                       @num.out @no=modssn2/@nmm={modmomin;modmomax}.
                       @sortv sortv=modssn2/sorttp=1.
                   end if.
                   compute modssn2n=nrow(modssn2).
                   compute modxssv={make(modssn2n,1,1),
                        make(modssn2n,1,(2*modtm)),make(modssn2n,1,0),modssn2}.
                   do if (modm=8).
                       compute modxssv={modxssv,modxssv(:,2)&*modssn2}.
                   end if.
                   compute modssn3=make(modssn2n,8,0).
                   compute modssn3(:,1)=make(modssn2n,1,modtm).
                   compute modssn3(:,2)=modssn2.
                   compute modssn3(:,3)=modxssv*modb.
                   compute modssn3(:,4)=sqrt(diag(modxssv*modbcov*t(modxssv))).
                   compute modssn3(:,5)=modssn3(:,3)&/modssn3(:,4).
                   do if (dichy=2).
                       compute modssn3(:,6)=2*(1-cdfnorm(abs(modssn3(:,5)))).
                       else.
                       compute modssn3(:,6)=2*(1-tcdf(abs(modssn3(:,5)),rg1df2)).
                   end if.
                   compute modssn3(:,7:8)={(modssn3(:,3)-modssn3(:,4)*qtout),
                        (modssn3(:,3)+modssn3(:,4)*qtout)}.
                   print modssn3/title '====================> Johnson-Neyman'
                          /cname=modssnm/rname=modssrnm/format=!dcms.
     do if (!modjnp=1).
         print /title '====================> Code for Plot (Johnson-Neyman)'.
         print /title=!plonote/format=a8/space=0.
         print modjnptt/title='  '/format=a8/space=0.
         print {modssn3(:,1:3),modssn3(:,7:8)}/title 'Begin data.'/format=!dcms/space=0.
         print /title 'End data.'/format=a8/space=0.
         print modjnpt2/title="Graph/scatterplot(overlay)="/format=a8/space=0.
     end if.
               end loop.
            end if.
        else if (modm=9 or modm=10).
            compute modssn2=modssn.
            compute modb=rgb(2:5,1).
            compute modbcov=rgbcov(2:5,2:5).
            do if (modm=10).
                compute modb=rgb(2:6,1).
                compute modbcov=rgbcov(2:6,2:6).
            end if.
            do if (modm=9).
                compute modjna=qtout2*modbcov(4,4)-modb(4,1)**2.
                compute modjnb=2*qtout2*modbcov(1,4)-2*modb(1,1)*modb(4,1).
                compute modjnc=qtout2*modbcov(1,1)-modb(1,1)**2.
                @eq.slv @eq.num={modjna,modjnb,modjnc}.
                do if (@eq.rtc=1).
                    compute modssn2={modssn2;@eq.rt}.
                    @num.out @no=modssn2/@nmm={modmomin;modmomax}.
                    @sortv sortv=modssn2/sorttp=1.
                end if.
                else.
                compute modjna=qtout2*modbcov(5,5)-modb(5,1)**2.
                compute modjnb=2*qtout2*modbcov(4,5)-2*modb(4,1)*modb(5,1).
                compute modjnc=2*qtout2*modbcov(1,5)+qtout2*modbcov(4,4)-
                    2*modb(1,1)*modb(5,1)-modb(4,1)**2.
                compute modjnd=2*qtout2*modbcov(1,4)-2*modb(1,1)*modb(4,1).
                compute modjne=qtout2*modbcov(1,1)-modb(1,1)**2.
                @eq.slv @eq.num={modjna,modjnb,modjnc,modjnd,modjne}.
                do if (@eq.rtc=1).
                    compute modssn2={modssn2;@eq.rt}.
                    @num.out @no=modssn2/@nmm={modmomin;modmomax}.
                    @sortv sortv=modssn2/sorttp=1.
                end if.
            end if.
            compute modssn2n=nrow(modssn2).
            compute modssn2={modssn2;modpapv}.
            compute modssn2m=nrow(modssn2).
            compute modxssv={make(modssn2m,1,1),make(modssn2m,2,0),modssn2}.
            compute modpapvn=nrow(modpapv).
            do if (modm=10).
                compute modxssv={modxssv,modxssv(:,4)&**2}.
                do if (modpapvn>1).
                    loop #m10=1 to modpapvn.
                        loop #m11=1 to modpapvn.
                            do if (#m10<#m11).
                                compute modxssv={modxssv;modxssv(modssn2n+#m10,:)-
                                    modxssv(modssn2n+#m11,:)}.
                            end if.
                        end loop.
                    end loop.
                end if.
            end if.
            compute modssn2l=nrow(modxssv).
            compute modssn=make(modssn2l,7,-99).
            compute modssn(:,1)=modxssv(:,4).
            compute modssn(:,2)=modxssv*modb.
            compute modssn(:,3)=sqrt(diag(modxssv*modbcov*t(modxssv))).
            compute modssn(:,4)=modssn(:,2)&/modssn(:,3).
            do if (dichy=2).
                compute modssn(:,5)=2*(1-cdfnorm(abs(modssn(:,4)))).
                else.
                compute modssn(:,5)=2*(1-tcdf(abs(modssn(:,4)),rg1df2)).
            end if.
            compute modssn(:,6:7)={(modssn(:,2)-modssn(:,3)*qtout),
                (modssn(:,2)+modssn(:,3)*qtout)}.
            compute modssnm2={nmx2(1,1),modssnm}.
            print modssn((modssn2n+1):modssn2m,:)/title '====================> Pick A Point '
                    /cname=modssnm/rname=modssrnm/foramt=!dcms.
            do if (modpapvn>1 and modm=10).
                compute tdata1=make(modpapvn*(modpapvn-1)/2,2,0).
                compute tdata2=0.
                loop #md1= 1 to modpapvn.
                    loop #md2=1 to modpapvn.
                        do if (#md1<#md2).
                            compute tdata2=tdata2+1.
                            compute tdata1(tdata2,:)={modssn(modssn2n+#md1,2),
                            modssn(modssn2n+#md2,2)}.
                        end if.
                    end loop.
                end loop.
                compute tdata1={tdata1,modssn((modssn2m+1):modssn2l,2:7)}.
                print tdata1/title '====================> Slope Difference Test'
                        /cname=modssnm3/rname=modssrnm/format=!dcms.
            end if.
            do if (!modjn=1).
                print modssn(1:modssn2n,:)/title '====================> Johnson-Neyman'
                        /cname=modssnm/rname=modssrnm/format=!dcms.
            end if.
  do if (!modjnp=1).
      print /title '====================> Code for Plot (Johnson-Neyman)'.
      compute modjnptt={'Data ','list ','free/',nmx2(1,2),' Slope ','LLCI','ULCI.'}.
      print /title=!plonote/format=a8/space=0.
      print modjnptt/title='  '/format=a8/space=0.
      print {modssn(1:modssn2n,1:2),modssn(1:modssn2n,6:7)}
          /title 'Begin data.'/format=!dcms/space=0.
      print /title 'End data.'/format=a8/space=0.
      compute tdata1={nmx2(1,2),nmx2(1,2),nmx2(1,2),'With';'Slope',' LLCI','ULCI','(Pair).'}.
      print tdata1/title="Graph/scatterplot(overlay)="/format=a8/space=0.
  end if.
        end if.
        do if (!modsave=1).
            print /title '*************************************** Save Data '+
            '***************************************'.
            save {eregid,edty,modx}/outfile=*/vars=EReg_ID !e_y .
            compute tdata1=ncol(modx)+1.
            compute prnm=t(colnm(1,3:(tdata1+1))).
            compute prnm={prnm,t(modxnm(1,2:tdata1))}.
            print prnm/title'>>>>> Column Names vs. Variables' 
                /cname={'COL Name','Variable'}  /rname={'List   :'}/format=a8.
        end if.
    end if.
end if.
do if (emod=5).
    do if (dichy=2).
        compute ero={ero;15}.
    end if.
    do if (sero=0).
        do if (nvarx2>(ncase-2)).
            compute ero={ero;7}.
        end if.
    end if.
    compute rrbtn=!rrbtn.
    do if (rrbtn>0).
        print rrbtn/title'Bootstrap'/rname={'Samples:'}/format=f10.0.
        do if (!rrbtm=1).
        do if (!quote(!seed)="random").
        print {'Random'}/title='Type   : Percentile'/rlabel='Seed   :'/format=a8/space=0.
        else.
        print seed/title='Type   : Percentile'/rlabel='Seed   :'/format=a8/space=0.
        end if.
        else.
        do if (!quote(!seed)="random").
        print {'Random'}/title='Type   : Bias-corrected'/rlabel='Seed   :'/format=a8/space=0.
        else.
        print seed/title='Type   : Bias-corrected'/rlabel='Seed   :'/format=a8/space=0.
        end if.
        end if.
    end if.
    compute sero=csum(ero).
    do if (sero=0). 
        compute rrkmin=!rrkmin.
        compute rrkmax=!rrkmax.
        compute rrkinc=!rrkinc.
        do if(rrkmin<0). 
            compute rrkmin=0.
            compute rrkcrt=1.
        end if .
        do if (rrkmax<=rrkmin).
            compute rrkmin=0.
            compute rrkmax=1.
            compute rrkcrt=1.
        end if.
        compute tdata1=rrkmax-rrkmin.
        do if (rrkinc>tdata1 or rrkinc<=0).
            compute rrkinc=tdata1/10.
            compute rrkcrt=1.
        end if.
        compute rrkiter=trunc(tdata1/rrkinc)+1.
        compute tdata1=mod(tdata1,rrkinc).
        do if (tdata1>0).
            compute rrkiter=rrkiter+1.
        end if.
        loop #rrki=1 to rrkiter.
            do if (#rrki=1).
                compute rrks=rrkmin.
                else if (#rrki=rrkiter).
                compute rrks={rrks,rrkmax}.
                else.
                compute rrks={rrks,(rrkmin+(#rrki-1)*rrkinc)}.
            end if.
        end loop.
        compute tdata1={edty,edtx2}.
        @RR3 @RRyx=tdata1/@RRk=rrks/@RRext=1/@RRkest=1.
        compute tdata3={nmy,nmx2,'       .'}.
        do if (!rrdes=1).
            compute rrmin=cmin(tdata1).
            compute rrmax=cmax(tdata1).
            compute tdata2={'Mean', 'Std.Dev', 'Min', 'Max'}.
            print {t(@RRyxm),@RRsd,t(rrmin),t(rrmax)}/rname=tdata3/format=!dcms
                /cname=tdata2/title='******************************'+
                    '*** Descriptive Statistics ********************************'.
        end if.
        do if (!rrcor=1).
            print @RRyxr/title='********************************** Correlation Matrix ******'+
                '*****************************' /cname=tdata3/rname=tdata3/format=!dcms.
        end if.
        do if (!rrkfit=1).
            compute tdata2={'K','R-square','adj.R-sq','MSE','AIC','BIC'}.
            print {t(rrks),t(@RRmfit)}/title='************************************ '+
                    'K vs. Model Fit ************************************'
                    /cname=tdata2/rname={'Results:'}/format=!dcms.
        end if.
        compute tdata3={'K',nmx2}.
        do if (!rrkstd=1).
            print {t(rrks),t(@RRbstd)}/cname=tdata3/rname={'Results:'}/format=!dcms
                /title='******************************** K vs. Std.Coefficients '+
                    '*********************************'.
        end if.
        do if (!rrkvif=1).
            print {t(rrks),t(@RRvif)}/cname=tdata3/rname={'Results:'}/format=!dcms
                /title='*************************************** K vs. VIF '+
                    '***************************************'.
        end if.
        compute tdata4=!rrplot1+!rrplot2+!rrplot3.
        do if(tdata4>0) .
            print /title=!plotitle.
            print /title=!plonote/space=0.
            do if (!rrplot1=1).
                print /title'====================> K vs. Model Fit'/format=a8.
                print /title 'Data list free/K Rsquare AdjRsq MSE AIC BIC .'/format=a8.
                print {t(rrks),t(@RRmfit)}/title 'Begin data.'/format=!dcms/space=0.
                print /title 'End data.'/format=a8/space=0.
                print /title 'Graph/scatterplot=K with Rsquare.'/foramt=a8/space=0.
                print /title 'Graph/scatterplot=K with MSE.'/foramt=a8/space=0.
                print /title 'Graph/scatterplot(overlay)=K K with AIC BIC(pair).'/space=0.
            end if.
                !let !tt = !concat('K')
                !let !tt2 = !concat('Graph/line(multiple)=Value','(')
                !do !vars !in (!e_x)
                    !let !tt = !concat(!tt,!blanks(1),!vars)
                    !let !tt2 = !concat(!tt2,!blanks(1),!vars)
                !doend 
                !let !tt = !concat('Data list free/',!tt)
                !let !tt = !concat(!tt,' .')
                !let !tt =!quote(!tt)
                !let !tt2 = !concat(!tt2,' ) BY K.')
                !let !tt2 =!quote(!tt2)
            do if (!rrplot2=1).
                print /title'====================> K vs. Std.Coefficients'/format=a8.
                print /title=!tt /foramt=a8.
                print {t(rrks),t(@RRbstd)}/title 'Begin data.'/format=!dcms/space=0.    
                print /title 'End data.'/format=a8/space=0.
                print /title=!tt2/format=a8/space=0.
            end if.
            do if (!rrplot3=1).
                print /title'====================> K vs. VIF'/format=a8.
                print /title=!tt /foramt=a8.
                print {t(rrks),t(@RRvif)} /title 'Begin data.'/format=!dcms/space=0.    
                print /title 'End data.'/format=a8/space=0.
                print /title=!tt2/format=a8/space=0.
            end if.
        end if.
        print @RRks/title='****************************** K from Different Estimators ***'+
            '***************************'/cname=@RRkscnm/rname=@RRksrnm/format=!dcms.
        compute rrkm=rnd(!rrkm).
        do if (rrkm<1 or rrkm>27).
            compute rrkfin=!rrkfin(1,1).
            else.
            compute rrkfin=@RRks(rrkm,1).
        end if.
        @RR3 @RRyx=tdata1/@RRk=rrkfin/@RRext=1/@RRkest=0.
        compute rrbrnm={'Constant';t(nmx2)}.
        compute rrbcnm1={'Coeff','S.E.','t','p','Std.Coef','VIF'}.
        print rrkfin/title'************************************ Ridge Rgression '+
            '***********************************'/rname={'K value:'}/format=!dcms.
        print {t(@RRmtest)}/title '>>>>> Ridge Model Test'
            /rname={'Results:'}/cname={'F.v','F.df1','F.df2','F.p','Coeff.df'}/format=!dcms. 
        print {t(@RRmfit)}/title '>>>>> Ridge Model Fit'
            /rname={'Results:'}/cname={'R-square','Adj.R-sq','MSE','AIC','BIC'}/format=!dcms. 
        print {@RRbr,@RRbrse,@RRbt,@RRbp}/title= '>>>>> Ridge '+
            'Coefficients (scaling predictors in correlation form)'
            /cname=rrbcnm1/rname=rrbrnm/format=!dcms.
        compute xx={make(ncase,1,1),edtx2}.
        compute tpred=xx*@RRb.
        compute tyss=cssq(edty-csum(edty)/ncase).
        compute tess=cssq(edty-tpred).
        compute trss=tyss-tess.
        compute trsq=trss/tyss.
        compute tdata2={0,0;@RRbstd,@RRvif}.
        print {@RRb,@RRbse,@RRbt,@RRbp,tdata2}/title= '>>>>> Coefficients'
            /cname=rrbcnm1/rname=rrbrnm/format=!dcms.
        compute rrbraw=@RRb.
        do if (rrbtn>0).
            compute rrbbt=make(rrbtn,(nvarx2+1),-999).
            compute btnero=1.
            compute rrbtdt=tdata1.
            loop #rbtn=1 to rrbtn.
                loop #rbtn2=1 to ncase.
                    compute rrbtdt(#rbtn2,:)=tdata1(trunc(uniform(1,1)*ncase+1),:).
                end loop.
                @candr crdata=rrbtdt.
                do if (candrn>0).
                    compute btero=btero+1.
                    else.
                    @RR3 @RRyx=rrbtdt/@RRk=rrkfin/@RRext=0/@RRkest=0.
                    compute rrbbt(btnero,:)=t(@RRb).
                    compute btnero=btnero+1.
                end if.
            end loop.
            do if (btero>0).
                loop #rbtn3=1 to btero.
                    compute tdata2=trunc(uniform(1,1)*btnero)+1.
                    compute rrbbt((rrbtn+1-#rbtn3),:)=rrbbt(tdata2,:).
                end loop.
            end if.
            loop #rbtn4=1 to (nvarx2+1).
                @sortv sortv=rrbbt(:,#rbtn4).
            end loop.
            @btsot btrawd=t(rrbraw)/btalld=rrbbt/btcitp=!rrbtm/btnewnm=xdm1
                /btcih=eregcih/btcil=eregcil.
            compute xdm1=t(xdm1).
            compute tdata2={rrbraw,xdm1(:,2),xdm1(:,4:6)}.
            compute tdata1={'Coeff','BootMean','BootSE','BootLLCI','BootULCI'}.
            print tdata2/title= '>>>>> Coefficients(Bootstrap)'/rname=rrbrnm
                /format=!dcms/cname=tdata1.
        end if.
    end if.
end if.
/** check errors **/.
compute sero=csum(ero)+btero+rrkcrt.
do if (sero>0).
    print /title='************************************ Errors and Notes '+
    '***********************************'/format=a8.
    do if (rrkcrt=1).
        print /title='Note   : Unreasonable K values had been corrected.'/format=a8.
    end if.
    do if (btero>0).
        print btero/title'Note   : Due to estimation problems, some bootstrap samples had to be '+
    'replaced.' /rname={'Size   :'}/format=!dcms.
    end if.
    compute nero=nrow(ero).
    loop #r=1 to nero.
        do if (ero(#r,1)=1).
            print /title='Error  : Some variables are constant.'/format=a8 .
            else if (ero(#r,1)=2).
            print /title='Error  : Two variables are the same.'/format=a8 .
            else if (ero(#r,1)=3).
            print /title='Error  : Singular or near singular data matrix.'/format=a8 .
            else if (ero(#r,1)=5).
            print /title='Error  : The number of the IVs must be greater than 1 in principal '+
    'component regression.'/format=a8 .
            else if (ero(#r,1)=12).
            print /title='Note   : Unreasonable confidence intervals have been '+
    'corrected.'/format=a8.
            else if(ero(#r,1)=7).
            print /title'Error  : The sample size is too small.'/format=a8.
            else if(ero(#r,1)=6).
            print /title'Error  : The number of IVs cannot exceed 100 in PCR.'/format=a8.
            else if(ero(#r,1)=10).
            print /title'Error  : The number of knots must be between 1 and 99.'/format=a8.
            else if(ero(#r,1)=11).
            print /title'Error  : The Q values are out of the data range.'/format=a8.
            else if(ero(#r,1)=13).
            print /title'Error  : The number of knots must be between 3 and 99'+
               ' in restricted cubic spline model.'/format=a8.
            else if(ero(#r,1)=14).
            print /title'Error  : You have not specified a moderator in this model.'/format=a8.
            else if(ero(#r,1)=15).
            print /title'Error  : This analysis does not allow dichotomous DV.'/format=a8.
        end if.
    end loop.
end if.
END MATRIX.
!enddefine. 
ereg1.0.1 emod=%%module%%
/e_y=  %%y%%
/e_x=%%x%%
/ci=%%ci%%
/seed=%%seed%%
/rgpm=%%rgpm%%
/rgpdn=%%rgpdn%%
/rgprbst=%%rgprbst%%
/rgpmol=%%rgpmol%%
/rgpmov=%%rgpmov%%
/rgpemp=%%rgpemp%%
/rgpemq=%%rgpemq%%
/mrres=%%mrres%%
/mrdsts=%%mrdsts%%
/mrdcasen=%%mrdcasen%%
/mrcoefev=%%mrcoefev%%
/mrca=%%mrca%%
/mrda=%%mrda%%
/mrallsub=%%mrallsub%%
/mrassot=%%mrassot%%
/mrdes=%%mrdes%%
/mrcovb=%%mrcovb%%
/mrheter=%%mrheter%%
/mrrbste=%%mrrbste%%
/mrplot1=%%mrplot1%%
/mrplot2=%%mrplot2%%
/qedttp=%%qedttp%%
/qerbst=%%qerbst%%
/qepap=%%qepap%%
/qepapv=%%qepapv%%
/qejn=%%qejn%%
/qejnp=%%qejnp%%
/qesave=%%qesave%%
/rrkmin=%%rrkmin%%
/rrkmax=%%rrkmax%%
/rrkinc=%%rrkinc%%
/rrkstd=%%rrkstd%%
/rrkfit=%%rrkfit%%
/rrkvif=%%rrkvif%%
/rrkm=%%rrkm%%
/rrkfin=%%rrkfin%%
/rrbtn=%%rrbtn%%
/rrbtm=%%rrbtm%%
/rrplot1=%%rrplot1%%
/rrplot2=%%rrplot2%%
/rrplot3=%%rrplot3%%
/rrdes=%%rrdes%%
/rrcor=%%rrcor%%
/pcrdes=%%pcrdes%%
/pcrcor=%%pcrcor%%
/pcrallrb=%%pcrallrb%%
/pcreva=%%pcreva%%
/pcreve=%%pcreve%%
/pcrcm=%%pcrcm%%
/pcrcsm=%%pcrcsm%%
/pcrplot1=%%pcrplot1%%
/pcrplot2=%%pcrplot2%%
/pcrplot3=%%pcrplot3%%
/pcrtype=%%pcrtype%%
/pcrnum=%%pcrnum%%
/pcrbtn=%%pcrbtn%%
/pcrbtm=%%pcrbtm%%
/spltp=%%spltp%%
/splktp=%%splktp%%
/splqv=%%splqv%%
/splrbst=%%splrbst%%
/splsave=%%splsave%%
/splplot=%%splplot%%
/splplotn=%%splplotn%%
/splemmp=%%splemmp%%
/modm=%%modm%%
/moddttp=%%moddttp%%
/modsave=%%modsave%%
/modrbst=%%modrbst%%
/modpap=%%modpap%%
/modpapv=%%modpapv%%
/modjn=%%modjn%%
/modjnp=%%modjnp%%
/modjn2=%%modjn2%%
/modjnx=%%modjnx%%
/dcms=%%dcms%%
.
RESTORE.
