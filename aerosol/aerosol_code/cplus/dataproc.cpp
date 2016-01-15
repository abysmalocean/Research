#include "dataproc.h"

int dataproc::readCaseCnt(string line)
{
	if(line.empty())
		return 1;
	else
		return  atoi(line.c_str());
	return 0;
}
    
void dataproc::readCase(string line)
{	
	        
	if(line.empty())
    	return;
	int linecnt = 0;
	int csnum=0;
	string strcsnm;
	while(!fin.eof())
	{	
		getline(fin,line);
		linecnt++;
		if(line.empty()) 
			continue;
		getline(fin,line);
		strcsnm = line.substr(0,39);
		xa = atof(line.substr(40,10).c_str());
		xb = atof(line.substr(50,10).c_str());
		relerr = atof(line.substr(60,10).c_str());
		abserr = atof(line.substr(70,10).c_str());
		getline(fin,line);
		nvar   = atoi(line.substr(0,2).c_str());
		nmxvar = atoi(line.substr(2,2).c_str());
		ndiscr = atoi(line.substr(4,2).c_str());
		nmulgv = atoi(line.substr(6,2).c_str());
		nprnt  = atoi(line.substr(8,2).c_str());
		iswtch = atoi(line.substr(10,2).c_str());
		isstat = atoi(line.substr(12,2).c_str());
		igrid  = atoi(line.substr(20,3).c_str());
		nptbet = atoi(line.substr(30,1).c_str());
		ncsorc = atoi(line.substr(31,1).c_str());
		ncremv = atoi(line.substr(32,1).c_str());
		ndsorc = atoi(line.substr(33,1).c_str());
		ndremv = atoi(line.substr(34,1).c_str());
		icondn = atoi(line.substr(35,1).c_str());
		isplin = atoi(line.substr(36,1).c_str());
		isetup = atoi(line.substr(37,1).c_str());
		ispswt = atoi(line.substr(38,1).c_str());
		isetsw = atoi(line.substr(39,1).c_str());
		itimdp = atoi(line.substr(40,1).c_str());
		for(int j=0; j<10; j++ )
		{
			iparm[j] = atoi(line.substr(50+j,1).c_str()); 
		}
							
		issflg = 0;
		if( !(ndiscr>0) )
		{
			nmulgv = 0;
			ndsorc = 0;
			ndremv = 0;		
			isstat = 0;
		}
		if(relerr<=0.0) relerr = 0.001;
		if(nvar<10) nvar=10;
		if(nmxvar<nvar) nmxvar=nvar;
		if(nprnt<0) nprnt=1;
		getline(fin,line);
		for( int i=0; i<nprnt; i++ )
		{
			tout[i] = atof(line.substr(i*10,10).c_str());
		}
		nvarm1 = nvar-1;
		nvarp1 = nvar+1;
		nvarp2 = nvar+2;
		nvart  = nvar+ndiscr;
		nmulp1 = nmulgv+1;
		nmulp2 = nmulgv+2;
		igdmp1 = abs(igrid)+nmulp1;
		igdmp2 = igdmp1+1;
		if( ndiscr != 0 )
		{
			ndcls[0] = ndiscr + 1;
			getline(fin,line);
			for(int i=1; i<nmulp1-1; i++)
				ndcls[i]=atof(line.substr(i*10,10).c_str());
			ndcls[nmulp1-1]=2*ndiscr;
			getline(fin,line);
			for(int j=nvarp1-1; j<nvart; j++)
				a[i]=atof(line.substr((i-nvarp1+1)*10,10).c_str());
			for( int j=0; j<ndiscr; j++ )
				x[nvar+j] = xa * j;
			for( int k=0; k<nmulp1; k++ )
				x[k]=xa*double(ndcls[k]);
			xa=x[1];	
		}
		xboxa = xb/xa;
		algxba = log10(xboxa);		
		x[0]=xa;
		x[nvar-1]=xb;
		w[0]=0.0;
		w[nvar-1]=1.0;
		nmulgv+=1;
		if(nmulgv!=0)
		{
			for(int i=1;i<nmulp1;i++)
				w[i]=log10(x[i]/xa)/algxba;
		}
	
		if(igrid == 0)
		{
			//igrid == 0
			for(int i=nmulp2; i<nvarm1; i++)
				w[i]=double(i-nmulp1)/double(nvar-nmulp1)*(1.0-w[nmulp1])+w[nmulp1];
		}
		else //igrid != 0, go to 16
		{  
			if(igrid>0)//igrid > 0
			{
				getline(fin,line);
				for(int i=nmulp2-1; i<igdpm1;i++)
				{
					x[i]=atof(line.substr((i-nmulp2+1)*10,10).c_str());
				}
			}
			if(igrid<0) //igrid < 0
			{
				getline(fin,line);
				for(int i=nmulp2-1; i<igdpm1; i++)
				{
					x[nvarm1+nmulp2-i-1]=atof(line.substr((i-numpl2+1)*10,10).c_str());
				}
			}
			int index;
			for(int j=nmulp2-1; j<igdpm1;j++)
			{
				index = j;
				if(igrid<0)
					index = nvarm1+nmulp2-j
				w[index] = log10(x[i]/xa)/algxba;
			}
			if(igdpm1>=nvarm1)
			{
			
			}
			if(igrid >0.0)
			{
				nvgdm1=nvar+igrid-1;
				for(int i=nmulp2-1; i<nvgdm1; i++)
				{
					w[i]=doubld(i-nmulp1)/double(nvar+igrid-nmulp1)*(w[nvar+igrid]-w[nmulp1])+w[nmulp1];
				}
			}
		}
			
		int index2=0;
		for(int i=igdpm2-1; i<nvarm1; i++)
		{
			index = i;
			if(igrid<0)
			{
				index = nmulp1+index-igdpm1;
				x[i]=xa*pow(xboxa,w[i]);
			}
		}
		//call xintl()
		for( int i=0; i<nvar; i++ )
			a[i]=a[i]*algxba*x[i];
		nquad2=nquadp/2;
		int index3=0;
		for( int i=0; i<nquadp; i++ )
		{
			index = max(nquad2-i+1,i-nquad2);
//			fa
/*				if(ndiscr!=0)
				{
					for( int i=nvarp1-1; i<nvart; i++ )
						a[i]=atof(line.substr((i-nvarp1+1)*10,10).c_str());
					for( int j=0; j<ndiscr; j++ )
						x[nvar+j] = xa * j;
					for( int k=0; k<nmulp1; k++ )
						x[k]=xa*double(ndcls[k]);
					xa=x[1];	
	    		}*/
		}
}
