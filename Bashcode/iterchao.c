#include<stdio.h>
FILE *fp;
main()
{
char s[200],t[200];
int i,j,k,c,XX,YY;
double xyz[385][301],xx,yy,zz,Xfac,Yfac;    /*fix the array size of slab grid*/
fp=fopen("slab1.0.grid","r");
for(i=300;i>-1;i--){
  for(j=384;j>-1;j--){
    fscanf(fp,"%lf %lf %lf",&xx,&yy,&zz);
    xyz[j][i] = zz;
    /*printf("%7.3lf %7.3lf %d %d %8.4lf\n",xx,yy,j,i,zz);*/
  }
}
fclose(fp);

fp=fopen("chat00p","r");
while(c!=EOF){
  for(k=0;k<6;k++){
    for(i=0;i<200&&(c=getc(fp))!=EOF&&c!='\n';i++) s[i]=c;
    if(c==EOF) continue;
    s[i] = '\0'; printf("%s\n",s);
  }

  for(i=0;i<200&&(c=getc(fp))!=EOF&&c!='\n';i++) s[i]=c;
  if(c==EOF) break;
  for(j=0;j<200;j++){
    for(i=0;i<200&&(c=getchar())!=EOF&&c!='\n';i++) s[i]=c;
    s[i] = '\0';
    if(s[1]=='0'&&s[3]=='-'&&s[6]=='-') break;
  }
  /*printf("%s\n",s);*/
  
  for(k=21;k<40;k++) if(s[k]==' ') s[k]='0';
  t[0]=s[1]; t[1]=s[2]; t[2]=s[4]; t[3]=s[5]; t[4]=s[7]; t[5]=s[8];
  t[6]=s[10]; t[7]=s[11]; t[8]=s[12]; t[9]=s[13]; t[10]=s[15]; t[11]=s[16];
  t[12]=s[18]; t[13]=s[19]; t[14]=s[20]; t[15]=s[21]; t[16]=s[22]; t[17]=' ';
  t[18]=s[24]; t[19]=s[25]; t[20]=s[27]; t[21]=s[28]; t[22]=s[30]; t[23]=s[31];
  t[24]=s[32]; t[25]=' '; t[26]=s[34]; t[27]=s[35]; t[28]=s[37]; t[29]=s[38];
  t[30] = '-'; t[31] = '\0';
  printf("%s",t);
  
  for(k=0;k<17;k++) if(t[k]==' ') t[k] = '0';
  xx = (t[22]-'0')*100.+(t[23]-'0')*10.+(t[24]-'0')*1.;     /*integer part of longitude, in deg */
  xx += ((t[26]-'0')*10.+(t[27]-'0')*1.+(t[28]-'0')*.1+(t[29]-'0')*.01)/60.;    /*decimal part of longitude, in deg*/
  yy = (t[15]-'0')*10.+(t[16]-'0')*1.;                      /*integer part of latitude, in deg */
  yy += ((t[18]-'0')*10.+(t[19]-'0')*1.+(t[20]-'0')*.1+(t[21]-'0')*.01)/60.;    /*decimal part of latitude, in deg*/
  XX = (int)(xx*100.)-12192;    /*(int)9.9==9, (int)9.3==9*/
  Xfac = xx*100.-12192.-(double)XX;     /*XX is the integer part, Xfac is the decimal part*/
  YY = (int)(yy*100.)-4650;
  Yfac = yy*100.-4650.-(double)YY;
  zz = xyz[(int)(xx*100.)-12192][(int)(yy*100.)-4650];

  zz = xyz[XX][YY]*(1.-Xfac)*(1.-Yfac) + xyz[XX][YY+1]*(1.-Xfac)*(Yfac)
     + xyz[XX+1][YY]*(Xfac)*(1.-Yfac) + xyz[XX+1][YY+1]*(Xfac)*(Yfac);

  printf("%4d\n",(int)(zz*100.));
  /*printf("%7.3lf %6.3lf   %d %d  %8.3lf  Xf= %6.3lf Yf = %6.3lf\n",xx,yy,(int)(xx*100.)-12150,(int)(yy*100.)-4650,zz,Xfac,Yfac);*/
}
fclose(fp);
}
