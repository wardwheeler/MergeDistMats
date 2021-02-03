{- |
Module      :  mergeDistMats.hs 
Description :  Progam to merge distance matrices with different leave sets for distance analysis 
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :  

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)


-}

module Main where

import System.IO
import System.Environment
import Data.List
import Data.Maybe
import qualified Data.Vector as V 
import Data.Matrix as M
import Text.ParserCombinators.Parsec
import Data.CSV
import Control.Parallel.Strategies
import Debug.Trace
import qualified Data.Number.Transfinite           as NT


-- | 
-- Map a function over a traversable structure in parallel
-- Preferred over parMap which is limited to lists
parmap :: Traversable t => Strategy b -> (a->b) -> t a -> t b
parmap strat f = withStrategy (parTraversable strat).fmap f


-- | getRawData  takes input ar of file name and returns rawData file contents
getRawData ::  [(String, Either ParseError [[String]])] -> [[[String]]]
getRawData inList =
    if null inList then []
    else 
      let (fileString, csvResult) = head inList
          rawData = case csvResult of
                    Left err -> error $ "Error parsing " ++ fileString ++ " " ++ show err 
                    Right result -> result
      in -- trace ("GRD:" ++ show rawData) 
      ((filter (/=[""]) rawData) : (getRawData $ tail inList))

-- | getMatrix takes tuple of leaf names and matrix checks for basic conditions (symmetry, square),
-- copnmverts to doubles and returnms list of [[Double]]
getMatrix :: ([String], [[String]]) -> M.Matrix Double
getMatrix (leafNames, rawData) =
  --trace ("in getMatrix:\n" ++ show leafNames ++ " " ++ (show $ tail rawData)) (
  if null leafNames then error "Empty leave name list"
  else if null rawData then error "Empty matrix data"
  else 
    if (length leafNames) /= (length $ tail rawData) then error ("Input matrix is not square: " ++ (show $ length leafNames) ++ " leaves and " ++ 
      (show $ length rawData) ++ " rows")
    else 
      let rowLengthCheck = foldl' (&&) True $ fmap (== length leafNames) $ fmap length (tail rawData)
      in
      if not rowLengthCheck then error "Row lengths do not equal leaf number"
      else
        let distMatrix = M.fromLists $ fmap (fmap (read :: String -> Double)) (tail rawData)
        in
        if not (M.isSymmetric distMatrix) then error ("Distance matrix is not symmetric--rows: " ++ (show $ rows distMatrix) ++ " columns: " ++ (show $ cols distMatrix))
        else distMatrix
        --)

-- | getSingleCell get average of cell bvalues for a pair of leaves through all matrices
getSingleCell :: String -> String -> [([String],M.Matrix Double)] -> Int -> Double -> Double
getSingleCell rowLeaf columnLeaf matrixList numMatrices curSum = 
  if null matrixList then 
      if numMatrices == 0 then trace ("Leaf pair " ++ rowLeaf ++ " " ++ columnLeaf ++ " never found in any matrices") (-1.0)
        -- error ("Leaf pair " ++ rowLeaf ++ " " ++ columnLeaf ++ " never found in any matrices")
      else curSum / (fromIntegral numMatrices)
  else 
      let (leafNames, distMatrix) = head matrixList
          rowIndex = elemIndex rowLeaf leafNames
          columnIndex = elemIndex columnLeaf leafNames
      in
      if rowIndex == Nothing then getSingleCell rowLeaf columnLeaf (tail matrixList) numMatrices curSum
      else if columnIndex == Nothing then getSingleCell rowLeaf columnLeaf (tail matrixList) numMatrices curSum
      else 
        let thisValue = distMatrix ! (fromJust rowIndex, fromJust columnIndex)
        in
        getSingleCell rowLeaf columnLeaf (tail matrixList) (numMatrices + 1) (curSum + thisValue)



-- | getCell determines a cell over all matrices taking average
getCells :: String -> [String] -> [([String],M.Matrix Double)] -> [Double]
getCells rowLeaf columnLeafList matrixList =
  if null columnLeafList then []
  else 
      let columnLeaf = head columnLeafList
          cellValue = getSingleCell rowLeaf columnLeaf matrixList 0 0.0 
      in
      cellValue : getCells rowLeaf (tail columnLeafList) matrixList

-- | makeRow take a list of rows left tot do (by leaf names) the full leaf list (for columns) 
-- and scans matrices for values taking average
makeRow :: [String] -> [String] -> [([String],M.Matrix Double)] -> [[Double]]
makeRow rowLeafList fullLeafList matrixList =
  if null fullLeafList then []
  else if null rowLeafList then []
  else 
      let rowLeaf = head rowLeafList
          thisRow = getCells rowLeaf fullLeafList matrixList
      in
      thisRow : makeRow (tail rowLeafList) fullLeafList matrixList

-- | mergeMatrices merges matrices by taking full leaf set and looks in each ,matrix in input list
-- for that pair's distance.  It also counts number of matrices a pair is in and takes the average of the 
-- distances.  Does the diagnoal and full square (not smart)
mergeMatrices :: [String] -> [([String],M.Matrix Double)] -> [[Double]]
mergeMatrices leafList matrixList =
  if null matrixList then error "No matrices to merge"
  else 
    makeRow leafList leafList matrixList

-- | normalizeMatrix takes a matrix and normaizes to [0,1]
normalizeMatrix :: M.Matrix Double -> M.Matrix Double
normalizeMatrix inMatrix =
  if inMatrix == M.empty then error "Empty matrix in normalizeMatrix"
  else
    let matrixRows = fmap (M.takeRow inMatrix) [0..((M.rows inMatrix) - 1)]
        maxVal = maximum $ fmap V.maximum matrixRows
    in
    if maxVal <= epsilon then inMatrix
    else M.map (/ maxVal) inMatrix

-- | epsilong for near zero comparison => 10^-9
epsilon :: Double
epsilon = (1.0 / 1000000000.0)

-- | nonZeroMinimum take a vector of doubles and returns smallest non-zero value
-- epsilon to  zero-
nonZeroMinimum :: [Double] -> Double -> Double 
nonZeroMinimum inList currentMinimum =
  if null inList then 
    if currentMinimum == NT.infinity then error ("Min value is infinity")
    else currentMinimum
  else 
    let firstVal = head inList
    in
    -- trace (show firstVal) (
    if NT.isNaN firstVal then nonZeroMinimum (tail inList) currentMinimum
    else if firstVal > currentMinimum then nonZeroMinimum (tail inList) currentMinimum
    else if firstVal <= epsilon then nonZeroMinimum (tail inList) currentMinimum
    else nonZeroMinimum (tail inList) firstVal
    --)

-- | integerizeMatrix takes input matrix and converts to integers by dividing all
-- elements with smallest non-zero value and rounding the result
integerizeMatrix ::  Double -> [[Double]] -> [[Int]]
integerizeMatrix precision inMatrix = 
  if null inMatrix then error "Empty matrix in integerizeMatrix"
    else
      let flatMat = concat inMatrix
          minNonZero = (nonZeroMinimum flatMat NT.infinity) / precision
      in
      --- trace ("Minimum is: " ++ show minNonZero  ++ " " ++ show precision ++ " " ++ show (nonZeroMinimum flatMat NT.infinity))
      fmap (fmap Prelude.round) $ fmap (fmap (/ minNonZero)) inMatrix

-- maximumDistance returns the sum of all the max distances in each matrix
maximumDistance :: Double -> [M.Matrix Double] -> Double
maximumDistance sumMaximum matrixList =
  if null matrixList then sumMaximum
  else 
    let thisMaximum' = V.maximum $ M.flatten $ head matrixList
        thisMaximum = if thisMaximum' <= epsilon then 1.0 else thisMaximum'
    in
    -- trace ("Maximum is: " ++ show thisMaximum) 
    maximumDistance (sumMaximum + thisMaximum) (tail matrixList)

-- | main driver
main :: IO ()
main = 
  do 
    -- Process arguments
    --  each csv file first line taxon/leaf names, subsequent lines are distances--must be symmetrical
    args <- getArgs
    if (length args < 2) then error ("Need first argument either 'normalize' (to scale matrix values on [0,1]) or 'weight' (to leave values as input]), and at least one csv files (usually > 1) to merge")
    else hPutStrLn stderr "Openning CSV files :"
    let normalize = if (args !! 0) == "normalize" then True else if (args !! 0) == "weight"  then False else error ("Parameter " ++ (args !! 0) ++ " not recongized")
    Prelude.mapM_ (hPutStrLn stderr) (fmap ("    " ++) args)

    -- Process input csv files
    -- rawDataLIst <- Prelude.mapM getRawData args
    rawFileList <- Prelude.mapM (parseFromFile csvFile) (tail args)
    let rawDataList = getRawData (Data.List.zip (tail args) rawFileList)
    -- hPutStrLn stderr ((show $ length rawDataList) ++ " " ) -- ++ show rawDataList)

    -- Check data integrity
    let leafNameList = Prelude.fmap head rawDataList 
    let matrixList = Prelude.fmap getMatrix (Data.List.zip leafNameList rawDataList) -- (fmap tail rawDataList))

    --Normalize matrices so [0,1] perhaps better for merging with same weights and missing data
    -- if not lets the input values accumulate
    let matrixListNorm = if normalize then Prelude.fmap normalizeMatrix matrixList else matrixList

    let fullLeafList = sort $ foldl' union (head leafNameList) (tail leafNameList)
    hPutStrLn stderr ("Complete leaf set: " ++ show fullLeafList)

    let mergedMatrix = mergeMatrices fullLeafList (Data.List.zip leafNameList matrixListNorm)

    -- Convert to integers by dividing each value by the non-zero minimum and rounding
    -- this to avoid floating point issues in Wag2020
    let maxDistance = maximumDistance 0.0 matrixList
    --hPutStrLn stderr ("Maximum overall distance: " ++ show maxDistance)
    let mergedMatrix' = integerizeMatrix maxDistance mergedMatrix

    let mergedMatsString = fmap (fmap show) mergedMatrix'

    -- Rename leaves if required
    -- Delete taxa if required

    hPutStrLn stdout (genCsvFile $ fullLeafList : mergedMatsString)

    hPutStrLn stderr "Done"


