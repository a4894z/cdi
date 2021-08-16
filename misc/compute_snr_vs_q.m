function [ snr_reconst, phot_px] = computeSNRvsQ( meas, sz, meas_region )


%#ok<*BDSCI>

%=========================================================================

%keep this for the time being, compare speeds for defining it here vs passing it
[ c, r ] = meshgrid( 1 : sz.c, 1 : sz.r );

%=========================================================================

%the q = 0 value:
jj = 1; kk = 1;

o_diam_circ = (((sz.c/2+1-c)./jj).^2 + ((sz.r/2+1-r)./jj).^2 < 1);

snr_reconst(kk) = sum(sum(o_diam_circ .* meas)) / (0E-9 + sum(sum(sqrt(o_diam_circ .* meas))));
phot_px(kk) = sum(sum(o_diam_circ .* meas));

%=========================================================================

jj_max = round(sz.c/2); 
annuls_smplng = 1;

%after q = 0, which corresponds to jj = 1, we now need to loop from:
loop_range = ((1 + annuls_smplng) : annuls_smplng : jj_max);

%=========================================================================

snr_reconst(2 : (1 + length(loop_range))) = 0;
phot_px(2 : (1 + length(loop_range))) = 0;

%{
snr_num(2 : (1 + length(loop_range))) = 0;
snr_denom(2 : (1 + length(loop_range))) = 0;
%}

%=========================================================================

%tic
kk = 2;
for jj = loop_range
  
  o_diam_circ = (((sz.c/2+1-c)./jj).^2 + ((sz.r/2+1-r)./jj).^2 < 1);
  i_diam_circ = (((sz.c/2+1-c)./(jj-1)).^2 + ((sz.r/2+1-r)./(jj-1)).^2 < 1);
  
  temp1 = meas_region .* (o_diam_circ - i_diam_circ);
  
  temp2 = temp1 .* meas;
  
  temp3 = sum(sum(temp2));
   
  phot_px(kk) = temp3 / sum(sum( temp1 ));
  
  %{
  snr_num(kk) = sum(sum(temp2));
  snr_denom(kk) = sum(sum(sqrt(temp2)));
  %}
  
  snr_reconst(kk) = temp3 / (eps + sum(sum(sqrt(temp2))));
  
  %{
  if isnan(snr_reconst(kk)),
    jj
  end
  %}
  
  kk = kk+1;
  
end
%toc

%snr_reconst = snr_num ./ snr_denom;
%snr_reconst(isnan(snr_reconst)) = 0.0707107;

%=========================================================================

%{
temp2 = 10*log10(snr_reconst); temp2 = temp2 .* (temp2>0);
figure; semilogx(temp2,'-or','LineWidth',2,'MarkerSize',4);  
xlim([12 length(snr_reconst)])
ylim([0 40])
%}

%=========================================================================

end

