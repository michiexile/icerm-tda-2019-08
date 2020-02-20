%% Read in the rotating cube data file

load 'rotatingCube.dat';

%% Construct the stream

dd = DistanceData(rotatingCube);

%% You'll want to repeat this until you have an rmax below 5E6

lm = WitnessStream.makeRandomLandmarks(dd,50);
rmax = WitnessStream.estimateRmax(dd,lm)

%% Construct the witness stream

ws = Plex.WitnessStream(1000,2,rmax,lm,dd);

%% Compute homology and cohomology

tic; intH = Plex.computePersistentHomology(ws); toc
tic; intB = Plex.computePersistentHomologyBasis(ws); toc
tic; intC = Plex.computePersistentCohomology(ws); toc

%% Extract longlived cohomology class

classes = [];
for i=1:length(intC)
    if(intC(i).dimension == 1 && intC(i).end-intC(i).start > 1)
        classes = [classes intC(i)];
    end
end

%% Extract boundary, coboundary and vector

interval = classes(1);
classVector = Plex.chainVector(ws,interval.basis,2500000);
boundaryMatrix = Plex.boundaryMatrix(ws,2500000,1);
coboundaryMatrix = Plex.coboundaryMatrix(ws,2500000,2);