function [rate] = ADK_rate_new(adk0, efield)
% E_max = 9.5080e+09 als Entwicklungspunkt in der ADK-Rate
rate = adk0 .* abs(efield./9.5080e+09).^13;
end