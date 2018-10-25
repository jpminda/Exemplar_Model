# Categorization model: Fit exemtype model to data, using hill-climbing algorithm and RMSD (RMSE)
# Amanda Lien - Categorization Lab, Western University
# Last Revised: May 2017

# randomly chooses values that sum to a user-chosen value
# in our case, we chose n number of values (number of traits) to sum to 1
randomSum <- function( N, M, sd = 0.1 )
{
  vec <- rnorm(N, M/N, sd)
  vec / sum(vec) * M
}

# determines whether changing the sensitivity or weight parameter exceeds our defined limits
exceeds <- function ( increment, type, x ) # type = weight or sensitivity (0)
{
  error<-FALSE
  if (increment == 0) # decrease
  {
    if ( (type == 0 && x < 2) | (type > 0 && x < 0.01 ) )
    {error<-TRUE}
  }
  else # increment = 1 = increase
  {
    if ( (type == 0 && x > 14) | (type > 0 && x > 0.99 ) )
    {error<-TRUE}
  }
  return( error )
}

exemplars<-read.table("SampleDatasets/6dExem.txt") # text file containing exemplars
observations<-read.table("SampleDatasets/6dObs.txt") # text file containing observational data to be fit to
stim<-read.table("SampleDatasets/6dStm.txt") # text file containing stim

# determines number of subjects/trials, dimensions (traits)
numTraits<-dim( exemplars )[2]
numStm<-dim( stim )[1]
numObs<-dim( observations )[1]
numExemplars<-dim( exemplars )[1]

# assumes that the first half of training set is Category A and the rest is Category B
# also assumes that there is an equal amt of A and B trainer sets
exemA<-matrix( unlist(exemplars[ 1:(numExemplars/2), 1:numTraits]), numExemplars/2, numTraits )
exemB<-matrix( unlist(exemplars[ (numExemplars/2+1):numExemplars, 1:numTraits ] ), numExemplars/2, numTraits)
numExemA<-dim(exemA)[1]
numExemB<-dim(exemB)[1]

# insert experimental data into matrix
stmData<-matrix( unlist(stim[ 1:numStm, 1:numTraits]), numStm, numTraits )
obsData<-matrix( unlist(observations[ 1:numObs, 1:numStm ]), numObs, numStm )

# produces random inital weights
makeWeights<-TRUE
while (makeWeights)
{
  makeWeights<-FALSE
  weightTable<-matrix( randomSum( numTraits, 1.0 ), numObs , numTraits, byrow=TRUE )
  for ( row in seq(numObs) )
  {
    if ( length(which(weightTable[row, 1:numTraits] < 0)) != 0 )
    { makeWeights<-TRUE }
  }
}

# initializes list of RMSDs at value of 1000 and sensitivites to a random value (b/w 0-10)
RMSD<-array( 1000, dim = numObs )
sens<-array( sample( c(1:2), numObs, replace = TRUE), dim = numObs )
prob<-matrix( NA, numObs, numStm )
whatChange<-array( NA, dim = numObs ) # will later determine what change (weight/sensitivit) will be made to each trial

# the following determines and applies changes to be made to weight or sensitivity in order to minimize RMSD
# once counter reaches value (user should choose), the program stops trying to minimize RMSD
# each weight and sensitivity has an equal chance of being altered
counter<-0

while ( counter < 1000 ) # <-- user chooses when to stop
{
  newWeights<-weightTable
  newSens<-sens

  whatChange<-sample( 0:numTraits, numObs, replace = TRUE) # randomly determine what change to make for each trial

  for ( currentTrial in seq( numObs )) # for each trial/subject, apply sens/weight change
  {
    direction<-sample( c(0:1), 1 ) # determines whether to increase(1)/decrease(0)
    if ( whatChange[currentTrial] == 0 ) # if zero, change sensitivitity
    {
      if ( exceeds(direction, whatChange[currentTrial], sens[currentTrial]) )
      {break} # if the value exceeds the limits, move on to the next trial/subject
      else # otherwise, add or subtract 1 based on direction (0 = decrease, 1 = increase)
      { if (direction == 0)
        { newSens[currentTrial]<-sens[currentTrial] - 0.01 }
        else
        { newSens[currentTrial]<-sens[currentTrial] + 0.01 }
      }
    }
    # CHANGE WEIGHT:
    # if the change in weight surpasses limits, terminate the loop and move onto next trial
    # if weight change is possible, make change and randomly choose another weight to compensate with
    else
    {
      if ( !exceeds(direction, whatChange[currentTrial], newWeights[currentTrial, whatChange[currentTrial]]) )
      { if (direction == 0) { newWeights[ currentTrial, whatChange[currentTrial] ]<-newWeights[ currentTrial, whatChange[currentTrial] ] - 0.01 }
        else { newWeights[ currentTrial, whatChange[currentTrial] ]<-newWeights[ currentTrial, whatChange[currentTrial] ] + 0.01 }
        secondCounter<-0
        while ( secondCounter < 20 && secondCounter >= 0 ) # looks for and applies change to weight randomly chosen to compensate with
        {
          secondCounter<-secondCounter + 1
          changeTo<-sample(seq(numTraits)[!seq(numTraits) %in% whatChange[currentTrial]], 1) # which other weight will be used to compensate
          swapDirection<-seq(0,1)[!seq(0,1) %in% direction]
          if ( !exceeds(swapDirection, changeTo, newWeights[currentTrial,changeTo]) ) # ensures that the compensating weight change doesn't surpass limits
          {
            if (swapDirection == 1)
            {
              newWeights[ currentTrial, changeTo ]<-weightTable[ currentTrial, changeTo ] + 0.01
              secondCounter<-(-1)
            }
            else
            {
              newWeights[ currentTrial, changeTo ]<-weightTable[ currentTrial, changeTo ] - 0.01
              secondCounter<-(-1)
            }
          }
        }
        if ( secondCounter == 20 )
        {
          if ( direction == 0 )
          { newWeights[ currentTrial, whatChange[currentTrial] ]<-newWeights[ currentTrial, whatChange[currentTrial] ] + 0.01 }
          else
          { newWeights[ currentTrial, whatChange[currentTrial] ]<-newWeights[ currentTrial, whatChange[currentTrial] ] - 0.01 }
        }
      }
    }
  }

  testRMSD<-array( dim = numObs ) # will contain calculated RMSDs to compare to pre-existing RMSDs
  predictedProb<-matrix( NA, numObs, numStm) # will contain predicted probabilities of choosing category A based on model and new parameters

  # finds the weighted differences between all exemplars to experimental data and stores
  # by storage and calculations in multi-dimensional arrays, where each 3rd dimension represents a trainer
  distancesA<-array( dim = c(numStm, numTraits, numExemA) )
  distancesB<-array( dim = c(numStm, numTraits, numExemB) )
  for ( eachTrial in seq(numObs) )
  {
    for ( eachTrainer in seq(numExemA) ) # for each trainer
    {
      distancesA[ 1:numStm, 1:numTraits, eachTrainer ]<-abs( matrix( exemA[eachTrainer,1:numTraits], numStm, numTraits, byrow = TRUE ) - stmData ) * matrix(newWeights[ eachTrial, 1:numTraits ], numStm, numTraits, byrow = TRUE)
      distancesB[ 1:numStm, 1:numTraits, eachTrainer ]<-abs( matrix( exemB[eachTrainer,1:numTraits], numStm, numTraits, byrow = TRUE ) - stmData ) * matrix(newWeights[ eachTrial, 1:numTraits ], numStm, numTraits, byrow = TRUE)

      for ( row in seq(numStm) ) # calculates similarity of each stimuli to category A and B across all exemplars
      {
        distancesA[ row, 1, eachTrainer ]<-exp( -newSens[eachTrial]*sum(distancesA[row, 1:numTraits, eachTrainer]) )
        distancesB[ row, 1, eachTrainer ]<-exp( -newSens[eachTrial]*sum(distancesB[row, 1:numTraits, eachTrainer]) )
      }
    }

    for ( eachRow in seq(numStm) ) # calculates means of sums of all exemplars and stores in first column of first matrix
    {
      distancesA[ eachRow, 1, 1 ]<-sum(distancesA[ eachRow, 1, 1:numExemA ])
      distancesB[ eachRow, 1, 1 ]<-sum(distancesB[ eachRow, 1, 1:numExemB ])
    }

    # calculate and store probability of classifying stimulus as category A, and then RMSD
    predictedProb[ eachTrial, 1:numStm ]<-distancesA[ 1:numStm, 1, 1 ] / (distancesA[ 1:numStm, 1, 1 ] + distancesB[ 1:numStm, 1, 1 ])
    testRMSD[ eachTrial ]<-sqrt( mean((obsData[eachTrial, 1:numStm] - predictedProb[eachTrial, 1:numStm])^2) )
  }

  # if any newly calculated RMSD is less than what was previously stored, replace it
  # otherwise, add 1 to the counter
  if ( length(which(testRMSD < RMSD)) == 0 )
  { counter <- counter + 1 }
  else
  {
    for ( eachChange in which(testRMSD < RMSD) )
    {
      weightTable[ eachChange, 1:numTraits ]<-newWeights[ eachChange, 1:numTraits ]
      sens[ eachChange ] <- newSens[ eachChange ]
      prob[ eachChange, 1:numStm ]<-predictedProb[ eachChange, 1:numStm ]

    }
    RMSD[ which(testRMSD < RMSD) ] <- testRMSD[ which(testRMSD < RMSD) ]
    counter <- 0
  }
}

print(sens)
print(weightTable)
print(RMSD)
print(prob)

#write.table(sens,file = "ExemSens.txt")
#write.table(weightTable,file = "ExemWeights.txt")
#write.table(RMSD,file = "ExemRMSD.txt")
#write.table(prob,file = "ExemProbs.txt")
