function F=get_sampling_frequency(N,prompt,title,def)
if nargin<1,
    N=1;
end
if nargin<2,
    prompt=[{'Enter sampling frequency:'},...
            {'Enter sampling frequency for the output dignal:'}];
end
if nargin<3,
    title='Sampling Frequency';
end
if nargin<4,
    def= [{'500'},{'2000'}];
end
lines = 1;

n=1:N;
f=inputdlg(prompt(n),title,lines,def(n));
if isempty(f)
   F=[];
else 
   for i=1:N
      F(i)=str2num(f{i});
   end
end

