% This draws panel A of Appendix Figure 1 in Amin et al., 2020
figure
subplot(1,3,1)
hold on
line([0 6], [0 7]);
line([0 7], [0 6]);
axis equal, axis tight
xlim([0 7]);
ylim([0 7]);

subplot(1,3,2)
hold on
line([0 2], [0 1.5]);
line([0 1.5], [0 2]);
axis equal, axis tight
xlim([0 7]);
ylim([0 7]);

subplot(1,3,3)
hold on
line([0 1.25], [0 2.25]);
line([0 2.25], [0 1.25]);
axis equal, axis tight
xlim([0 7])
ylim([0 7])
