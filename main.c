#include <stdio.h>
#include <stdlib.h>

double sqrt(double);
double logBase(double,double);
double ln(double);
double pow(double, double);
double exp(double);
double absolute(double);
double sign(double);
double deriv(double(*)(double),double );
double integr(double (*)(double),double , double ,double );
double sin(double);
double cos(double);
double tan(double);
double sec(double);
double csc(double);
double cot(double);
double W(double);
int factorial(int);
double lnFactorial(int);
double round(double);
double limit(double (*)(double),double);
double func(double);
double mod(double, double);
double gamma(double);

double arcsin(double);

int main()
{
    char* str= "hello";
    str[0]='a';
    printf("%s",str);
    return 0;
}


/*
Square root algorithm uses a binary search
*/
double sqrt(double x){
    double guess=x/2;
    double max=x;
    if(x<1){
        max=1;
    }
    double min=0;


    double prevGuess=guess;
    while(guess*guess!=x){
        prevGuess=guess;
        double guess2=guess*guess;
        if(guess2>x){
            max=guess;
            guess=(guess+min)/2;
        }
        else if(guess2<x){
            min=guess;
            guess=(guess+max)/2;
        }
        if(guess==prevGuess){
            break;
        }
    }
    return guess;

}



/*
uses the sqrt function to repeatedly divide
*/
double logBase(double base, double x){
    double log=0;
    double scaleStep=base;
    double step = 1;

    if(x<1){
        return -logBase(base,1/x);
    }

    while(x!=1&&scaleStep>1){
        if(x/scaleStep>=1){
            x/=scaleStep;
            log+=step;
        }
        else{
            scaleStep=sqrt(scaleStep);
            step/=2;
        }

    }
    return log;
}

const double E = 2.718281828459045235360;
const double PI = 3.14159265;

double ln(double x){
    return logBase(E,x);
}



/*
uses the sqrt function to repeatedly multiply
*/
double pow(double base, double x){

    if(x<0){
        return 1.0/pow(base,-x);
    }



    double ans = 1;
    double step=1;
    while(x!=0){
        if(x-step>=0){
            x-=step;
            ans*=base;
        }
        else{
            step/=2;
            base=sqrt(base);
        }

    }

    return ans;
}


double exp(double x){
    return pow(E,x);
}




//returns absolute value
double absolute(double x){
if(x<0){
    return -x;
}
return x;
}



//sign function
double sign(double x){
if(x>0){
    return 1;
}
else if (x<0){
    return -1;
}
return 0;

}


//uses limit process to find derivative
double deriv(double (*function)(double),double x){
    double deltaX=1;
    double prevDeriv=0;
    double d=1.0;
    do {
        deltaX*=0.9;
        if(deltaX==0){
            break;
        }
        prevDeriv=d;


        d=((*function)(x+deltaX)-(*function)(x))/(deltaX);
    }
    while(d!=prevDeriv);
    return d;
}


//uses limit process to find integral
double integr(double (*function)(double),double start, double end,double deltaX){
    double inte=1;
    double prevInte=0;

        if(deltaX==0){
            return 0;
        }
        prevInte=inte;
        double sum = 0;
        for(double i = start; i<end; i+=deltaX){
            sum+=(*function)(i);
        }
        sum*=deltaX;
        inte=sum;

    return inte;

}





/*
Uses Sum and difference rules and half angle rule
Algotithm works similar to the log and pow functions
*/
double sin(double x){
    x=mod(x,2*PI);
    //printf("(%lf,%lf) \n",orig,x);
    double stepAngle=PI;
    double stepValue=0;
    double angle=0;
    double value=0;

    while(angle!=x){
        if(angle+stepAngle>x){

            double cos=sqrt(1-stepValue*stepValue)*sign(PI-mod(stepAngle+PI/2,2*PI));

            stepValue=sqrt((1-cos)/2);
            //printf("%lf \n",stepValue);
            stepAngle/=2;
        }
        else{


                double sinA=value;
                double cosA=sqrt(1-value*value)*sign(PI-mod(angle+PI/2,2*PI));
                double sinB=stepValue;
                double cosB=sqrt(1-stepValue*stepValue)*sign(PI-mod(stepAngle+PI/2,2*PI));

                angle+=stepAngle;
                value=sinA*cosB+cosA*sinB;
        }
    //printf("%lf,%lf \n",angle,value);

    }



    return value;


}

/*
Uses Sum and difference rules and half angle rule
Algotithm works similar to the log and pow functions
*/
double cos(double x){
    x=mod(x,2*PI);
    double angle=0;
    double value=1;
    double stepAngle=PI;
    double stepValue=-1;

    while(angle!=x){
            if(angle+stepAngle>x){
                stepAngle/=2;
                stepValue=sqrt((1+stepValue)/2)*sign(PI-mod(stepAngle+PI/2,2*PI));
            }
            else{
                double cosA=value;
                double cosB=stepValue;
                double sinA=sqrt(1-value*value)*sign(PI-mod(angle,2*PI));
                double sinB=sqrt(1-stepValue*stepValue)*sign(PI-mod(stepAngle,2*PI));

                value=cosA*cosB-sinA*sinB;
                angle+=stepAngle;

            }

    }
    return value;

}




double tan(double x){
    return sin(x)/cos(x);
}


double sec(double x){
    return 1/cos(x);
}

double csc(double x){
    return 1/sin(x);
}



double cot(double x){
    return cos(x)/sin(x);
}



/*
My attempt on making the Lambert W function (inverse of xln(x))
I am using newtons method (xn=x(n-1)-f(x(n-1)/f'(x(n-1))))
f(x)=xln(x)
f'(x)=(x)(1/x)+(ln(x))(1)
f'(x)=1+ln(x)
f'(x)=ln(x)+1
*/
double W(double x){
   double input = x;
    double guess = (x*ln(x)-input)/(ln(x)+1);



    double newGuess=-1;
    while(newGuess!=guess){
        guess=newGuess;
        newGuess=guess-(guess*ln(guess)-input)/(ln(guess)+1);
    }
    return guess;
}





int factorial(int x){
    int ans=1;
    for(int i=1; i<=x; i++){
        ans*=i;
    }
    return ans;
}



double lnFactorial(int x){
    double ans=0;
    for (int i=1; i<=x; i++){
        ans+=ln(i);
    }
    return ans;
}



double round(double x){
    if(x<0){
        x=x-0.5;
    }
    else{
        x=x+0.5;
    }
    return x-mod(x,1.0);
}



double limit(double (*function)(double),double x){
    double deltaX= 1;
    double error=0.000000000000001;
    while((*function)(x+deltaX)!=(*function)(x-deltaX)&&deltaX>error){
            double old = (*function)(x+deltaX);
            deltaX*=0.9;

            if(absolute(old-(*function)(x+deltaX))>absolute((*function)(x+deltaX))){
                if(sign((*function)(x+deltaX))!=(*function)(x-deltaX)){
                    printf("asymptote?");
                    return 0;
                }
                else if (sign((*function)(x+deltaX))==1){
                    printf("+ infinity?");
                    return 0;
                }
                else{
                    printf("- infinity?");
                    return 0;
                }
            }



            if(deltaX<=error&&(sign((*function)(x+0.0001))!=sign((*function)(x+0.0001))||absolute((*function)(x+0.000001)-(*function)(x-0.000001))>0.01)) {
                printf("osiclating?");
                return 0;
            }



    }
    return (*function)(x+deltaX);
}


//used to test the derivative, integral, and limit functions
double func(double x){
    if(x==3){
        return 5;
    }
    return x;
}

/*
my implimation of the mod function since C does not seem to allow its use on doubles/floats
I use the mod function to simplify the trig
*/
double mod(double num, double dem){
    if(dem==0){
        return 0;
    }
    double quot=num/dem;
    while(quot>1){
        quot-=1;
    }
    while(quot<-1){
        quot+=1;
    }
    if(quot==1){
        return 0;
    }
    return dem*quot;
}


double gammaInput=0;
double gammaI(double x){
    return pow(x,gammaInput-1)*exp(-x);
}
double gamma(double z){
    gammaInput=z;
    return integr(gammaI,0,100,0.001);
}





double arcsin(double x){
    double stepA=PI/2;
    double stepV=1;

    double angle=0;
    double value=0;

    if(x<0){
        stepA=-stepA;
        stepV=-stepV;
    }

    while(absolute(value-x)>0.0000001){
        double cosA=sqrt(1-stepV*stepV)*sign(PI-mod(stepA+PI/2,2*PI));
        double sinA=stepV;
        double sinB=value;
        double cosB=sqrt(1-value*value)*sign(PI-mod(angle+PI/2,2*PI));
        if(x>0&&sinA*cosB+sinB*cosA>x){
            stepV=sqrt((1-cosA)/2);
            stepA/=2;
        }
        else if(x<0&&sinA*cosB+sinB*cosA<x){
            stepV=sqrt((1-cosA)/2);
            stepA/=2;
        }
        else{
            value=sinA*cosB+sinB*cosA;
            angle+=stepA;
        }
        //printf("%lf\n",value);
    }
    return angle;
}




