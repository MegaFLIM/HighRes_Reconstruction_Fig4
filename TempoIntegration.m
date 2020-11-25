function T = TempoIntegration(b)
%performs temporal integration
T = squeeze(sum(b,3));
end
