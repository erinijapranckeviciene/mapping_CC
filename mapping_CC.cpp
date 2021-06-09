/* mapping_CC.cpp : Defines the entry point for the console application. */
/* UNIX version */

/* #include "stdafx.h" */
/* #include <math.h> */
#include "mapping_CC.h"

int main(int argc, char * argv[])
{
	FILE * file_pattern=NULL;
	FILE * file_sequence=NULL;
	FILE * file_output=NULL;
	char * input_pattern=NULL;
	char * input_sequence=NULL;
	char * output_map=NULL;
	char * line_sequence;
    char * line_temp;
	int window_smoothing=1;
	int sigma=0;
	int i=0;
	int j=0;
	int k=0;
	int l=0;
    int width_matrix=0;
    int length_matrix=0;
	int maximal_length=0;
	int length_sequence=0;
	int * matrix_type;
	double * matrix;
	double * gaussian_smoothing;
	double * sequence_matrix;
	double * sum_matrix;
	double * sum_sequence;

	/* load input parameters */

	for(i=1;i<argc;i++)
		if(argv[i][0]=='-'&&argv[i][1]=='m')
			input_pattern=argv[i+1];

	for(i=1;i<argc;i++)
		if(argv[i][0]=='-'&&argv[i][1]=='s')
			input_sequence=argv[i+1];

	for(i=1;i<argc;i++)
		if(argv[i][0]=='-'&&argv[i][1]=='o')
			output_map=argv[i+1];

	for(i=1;i<argc;i++)
		if(argv[i][0]=='-'&&argv[i][1]=='w')
			window_smoothing=char2int(argv[i+1]);

	for(i=1;i<argc;i++)
		if(argv[i][0]=='-'&&argv[i][1]=='e')
			sigma=char2int(argv[i+1]);

	if(input_pattern==NULL||input_sequence==NULL||window_smoothing<1||sigma<0)
	{
		print_help();
		return -1;
	}

/*  open input and output files */

	file_pattern=fopen(input_pattern,"r");
	if(file_pattern==NULL)
	{
		printf("cannot open file %s \n",input_pattern);
		return -1;
	}

	file_sequence=fopen(input_sequence,"r");
	if(file_sequence==NULL)
	{
		printf("cannot open file %s \n",input_sequence);
		return -1;
	}

	if(output_map!=NULL)
	{
		file_output=fopen(output_map,"w");
		if(file_output==NULL)
		{
			printf("cannot create file %s\n", output_map);
			return -1;
		}
	}

	/* read header of matrix */
	line_sequence=read_line(file_pattern);

	for(i=0;line_sequence[i]<=32&&line_sequence[i]!='\0';i++)
	  ;

	for(line_temp=line_sequence+i;line_temp[0]!='\0';line_temp=word_next(line_temp))
	  width_matrix++;

	matrix_type = (int*) malloc(sizeof(int)*width_matrix*3);

	j=0;
	for(line_temp=line_sequence+i;line_temp[0]!='\0';line_temp=word_next(line_temp))
	{
		for(k=0;line_temp[k]!='\0'&&line_temp[k]>32;k++)
			;
		matrix_type[j*3]=k;
		maximal_length=(j==0||maximal_length>k) ? k : maximal_length;
		matrix_type[j*3+1]=(line_temp[0]=='A'||line_temp[0]=='T'||line_temp[0]=='G'||line_temp[0]=='C'||
			line_temp[0]=='a'||line_temp[0]=='t'||line_temp[0]=='g'||line_temp[0]=='c') ? 0 :
			(line_temp[0]=='W'||line_temp[0]=='S'||line_temp[0]=='w'||line_temp[0]=='s') ? 1 :
			(line_temp[0]=='R'||line_temp[0]=='Y'||line_temp[0]=='r'||line_temp[0]=='y') ? 2 : 3;
		matrix_type[j*3+2]=hash(line_temp,matrix_type[j*3],matrix_type[j*3+1]);
		if(matrix_type[j*3+1]>2||matrix_type[j*3+1]<0)
		{
			printf("wrong header of matrix\n");
			return -1;
		}
	 	j++;
	}

	/* read matrix */
	length_matrix=0;
	while(!feof(file_pattern)&&line_sequence!=NULL)
	{
		length_matrix++;
		free(line_sequence);
		line_sequence=read_line(file_pattern);
	}
	length_matrix--;
	if(length_matrix<1)
	{
		printf("empty matrix of pattern\n");
		return -1;
	}
	rewind(file_pattern);
	matrix = (double*) malloc(sizeof(double)*width_matrix*length_matrix);

	line_sequence=read_line(file_pattern);
	line_sequence=read_line(file_pattern);
	j=0;
	while(!feof(file_pattern)&&line_sequence!=NULL)
	{
		for(k=0;line_sequence[k]!='\0'&&line_sequence[k]<=32;k++)
			;
		line_temp=line_sequence+k;
		for(i=0;i<width_matrix&&line_temp[0]!='\0';i++)
		{
			matrix[i*length_matrix+j]=char2double(line_temp);
			line_temp=word_next(line_temp);
		}
		free(line_sequence);
		line_sequence=read_line(file_pattern);
		j++;
	}
/*	for(i=0;i<length_matrix;i++)
	{
		for(j=0;j<width_matrix;j++)
			printf("%f ",matrix[i*width_matrix+j]);
		printf("\n");
	}
	return 0;*/


	sum_matrix = (double*) malloc(sizeof(double)*width_matrix*3);
	sum_sequence = (double*) malloc(sizeof(double)*width_matrix*4);
	for(i=0;i<width_matrix;i++)
		for(j=0;j<3;j++)
			sum_matrix[i*3+j]=0;
	for(i=0;i<width_matrix;i++)
		for(j=0;j<length_matrix;j++)
		{
			sum_matrix[i*3]+=matrix[i*length_matrix+j];
			sum_matrix[i*3+1]+=matrix[i*length_matrix+j]*matrix[i*length_matrix+j];
			sum_matrix[i*3+2]=sqrt((j+1)*sum_matrix[i*3+1]-sum_matrix[i*3]*sum_matrix[i*3]);
		}

	/* gaussian smoothing matrix */
	gaussian_smoothing = (double*) malloc(sizeof(double)*window_smoothing);
	for(i=0;i<window_smoothing;i++)
		gaussian_smoothing[i]=exp(double(-(((i-int(window_smoothing/2))*sigma*2/double(window_smoothing))*
		((i-int(window_smoothing/2))*sigma*2/double(window_smoothing))/2)));

	/* read sequence*/
	line_sequence=read_line(file_sequence);
	while(!feof(file_sequence)&&line_sequence!=NULL)
	{
	/* create sequence profile */
		for(i=0;line_sequence[i]!='\0';i++)
			;
		sequence_matrix = (double*) malloc(sizeof(double)*(i-maximal_length+1)*width_matrix);
		for(j=0;j<width_matrix;j++)
			for(k=0;k<i-matrix_type[j*3]+1;k++)
				sequence_matrix[(i-maximal_length+1)*j+k]=0;
		for(j=0;j<width_matrix;j++)
			for(k=0;k<i-matrix_type[j*3]+1;k++)
				if(matrix_type[j*3+2]==hash(line_sequence+k,matrix_type[j*3],matrix_type[j*3+1]))
					for(l=0;l<window_smoothing;l++)
						if(k-int(window_smoothing/2)+l>=0&&k-int(window_smoothing/2)+l<i-matrix_type[j*3]+1)
							sequence_matrix[(i-maximal_length+1)*j+k-int(window_smoothing/2)+l]+=gaussian_smoothing[l];

/*		for(j=0;j<width_matrix;j++)
			for(k=0;k<i-matrix_type[j*3]+1;k++)
				printf("%d %f\n",k, sequence_matrix[(i-matrix_type[j*3]+1)*j+k]);
		return 0;*/

		length_sequence=i;

		/* calculation CC: 1 - sum, 2 - sum^2, 3 - sigma, 4 - CC=x(i)*y(i) */
		for(k=0;k<length_sequence-length_matrix-maximal_length+1;k++)
		{
			if(output_map!=NULL)
				fprintf(file_output,"%d",k+int(length_matrix/2)+1);
			else
				printf("%d ",k+int(length_matrix/2)+1);
			for(i=0;i<width_matrix;i++)
				for(j=0;j<4;j++)
					sum_sequence[i*4+j]=0;
			for(i=0;i<width_matrix;i++)
			{
				for(j=0;j<length_matrix;j++)
				{
					sum_sequence[i*4]+=sequence_matrix[i*(length_sequence-maximal_length+1)+k+j];
					sum_sequence[i*4+1]+=sequence_matrix[i*(length_sequence-maximal_length+1)+k+j]*
						sequence_matrix[i*(length_sequence-maximal_length+1)+k+j];
					sum_sequence[i*4+3]+=sequence_matrix[i*(length_sequence-maximal_length+1)+k+j]*
						matrix[i*length_matrix+j];
				}
				sum_sequence[i*4+2]=sqrt(j*sum_sequence[i*4+1]-sum_sequence[i*4]*sum_sequence[i*4]);
				if(sum_matrix[i*3+2]>0&&sum_sequence[i*4+2]>0)
					sum_sequence[i*4+3]=(length_matrix*sum_sequence[i*4+3]-sum_matrix[i*3]*sum_sequence[i*4])/
						(sum_matrix[i*3+2]*sum_sequence[i*4+2]);
				else
					sum_sequence[i*4+3]=0;
				if(output_map!=NULL)
					fprintf(file_output," %7.4f",sum_sequence[i*4+3]);
				else
					printf(" %7.4f",sum_sequence[i*4+3]);
			}
			if(output_map!=NULL)
				fprintf(file_output,"\n");
			else
				printf("\n");
		}

		/* free sequence */
		free(sequence_matrix);
		free(line_sequence);
		line_sequence=read_line(file_sequence);
		j++;
	}

	/* free */
	free(matrix);
	free(sum_matrix);
	free(sum_sequence);
	fclose(file_pattern);
	fclose(file_sequence);
	if(output_map!=NULL)
		fclose(file_output);

	return 0;
}
/* ///////////////////////////// Printing help ///////////////////// */
void print_help()
{
	printf("------------------------------------------------------------\n");
	printf(" Mapping DNA sequence by pattern matrix CC                  \n");
	printf("------------------------------------------------------------\n");
	printf("input parameters:                                           \n");
	printf("-m - input pattern                                          \n");
	printf("-s - input sequence                                         \n");
	printf("-------------------------                                   \n");
	printf(" no mandatory parameters:                                   \n");
	printf("-o - output file                                            \n");
	printf("-w - window of gaussian smoothing (recomend odd 3 - 11)     \n");
	printf("-e - sigma (recomend 1 - 5, 0 average)                      \n");
	printf("BMI, Medical school, University of Ottawa, S.Hosid 7.07.2011\n");
}
/* /////////////////////////// read line from file ///////////// */
char * read_line(FILE * file_name)
{
	char * read_string;
	char * read_string_temp;
	int i;
	int length_string=1024;
	char check='\0';
	i=0;
	read_string= (char*) malloc(sizeof(char)*length_string);
	while((check>=32&&i>0||i==0)&&!feof(file_name))
	{
		check=fgetc(file_name);
		if(check>=32)
		{
			read_string[i]=check;
			i++;
			if(i==length_string)
			{
				length_string*=2;
				read_string_temp= (char*) malloc(sizeof(char)*length_string);
				for(i=0;i<length_string/2;i++)
					read_string_temp[i]=read_string[i];
				free(read_string);
				read_string=read_string_temp;
			}
		}
	}
	if(i==0)
	{
		free(read_string);
		return NULL;
	}
	read_string[i]='\0';
	return read_string;
}
/* //////////////////////// next word from string ////////////////// */
char * word_next(char * line)
{
	int i;
	for(i=0;line[i]!=' '&&line[i]!='\0'; i++) ;
	for( ; line[i]==' '&&line[i]!='\0'; i++) ;
	return line+i;
}
/* ////////////// hash function converts word to int /////////////// */
int hash(char * c, int l, int type)
{
	int LENGTH_ALPHABET=0;
	int i=0;
	int j=0;
	int k=0;
	if(type==0)
	{
		LENGTH_ALPHABET=4;
		for(j=0;j<l;j++)
		{
			if(c[j]=='a'||c[j]=='A')
				i=0;
			else if(c[j]=='t'||c[j]=='T')
				i=1;
			else if(c[j]=='g'||c[j]=='G')
				i=2;
			else if(c[j]=='c'||c[j]=='C')
				i=3;
			else
				return -1;
			k=k*LENGTH_ALPHABET+i;
		}
	}
	else if(type==1)
	{
		LENGTH_ALPHABET=2;
		for(j=0;j<l;j++)
		{
			if(c[j]=='a'||c[j]=='A'||c[j]=='t'||c[j]=='T'||c[j]=='w'||c[j]=='W')
				i=0;
			else if(c[j]=='g'||c[j]=='G'||c[j]=='c'||c[j]=='C'||c[j]=='s'||c[j]=='S')
				i=1;
			else
				return -1;
			k=k*LENGTH_ALPHABET+i;
		}
	}
	else if(type==2)
	{
		LENGTH_ALPHABET=2;
		for(j=0;j<l;j++)
		{
			if(c[j]=='a'||c[j]=='A'||c[j]=='g'||c[j]=='G'||c[j]=='r'||c[j]=='R')
				i=0;
			else if(c[j]=='t'||c[j]=='T'||c[j]=='c'||c[j]=='C'||c[j]=='y'||c[j]=='Y')
				i=1;
			else
				return -1;
			k=k*LENGTH_ALPHABET+i;
		}
	}
	return k;
}
/* //////////////////////// Integer from string //////////////////// */
int char2int(char * line)
{
	int i=0;
	int sum=0;
	int kof=1;
	int sign=0;
	if(line[0]=='-'||line[0]=='+')
		sign=1;
	for(i=sign; line[i]>47&&line[i]<58; i++) ;
	i--;
	sum+=(line[i]>47&&line[i]<58) ? line[i]-48 : 0;
	for(i--; i>=sign; i--)
	{
		kof*=10;
		sum+=(line[i]-48)*kof;
	}
	if(line[0]=='-')
		sum*=-1;
	return sum;
}
/* //////////////////////// Double from string ///////////////////// */
double char2double(char * line)
{
	int i=0;
	int j=0;
	double sum=0;
	int kof=1;
	int sign=0;
	if(line[0]=='-'||line[0]=='+')
		sign=1;
	for(j+=sign; line[j]>47&&line[j]<58; j++) ;
	i=j-1;
	sum+=(line[i]>47&&line[i]<58) ? line[i]-48 : 0;
	for(i--; i>=sign; i--)
	{
		kof*=10;
		sum+=(line[i]-48)*kof;
	}
	i=j+1;
	for(j++; line[j]>47&&line[j]<58; j++) ;
	kof=1;
	for(; i<j ; i++)
	{
		kof*=10;
		sum+=((float)line[i]-48)/kof;
	}
	if(line[0]=='-')
		sum*=-1.0;
	return sum;
}
