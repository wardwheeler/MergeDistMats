{-# LANGUAGE CPP #-}
{-# LANGUAGE NoRebindableSyntax #-}
{-# OPTIONS_GHC -fno-warn-missing-import-lists #-}
module Paths_mergeDistMats (
    version,
    getBinDir, getLibDir, getDynLibDir, getDataDir, getLibexecDir,
    getDataFileName, getSysconfDir
  ) where

import qualified Control.Exception as Exception
import Data.Version (Version(..))
import System.Environment (getEnv)
import Prelude

#if defined(VERSION_base)

#if MIN_VERSION_base(4,0,0)
catchIO :: IO a -> (Exception.IOException -> IO a) -> IO a
#else
catchIO :: IO a -> (Exception.Exception -> IO a) -> IO a
#endif

#else
catchIO :: IO a -> (Exception.IOException -> IO a) -> IO a
#endif
catchIO = Exception.catch

version :: Version
version = Version [0,1,0,0] []
bindir, libdir, dynlibdir, datadir, libexecdir, sysconfdir :: FilePath

bindir     = "/Users/ward/.cabal/bin"
libdir     = "/Users/ward/.cabal/lib/x86_64-osx-ghc-8.10.2/mergeDistMats-0.1.0.0-inplace-mergeDistMats"
dynlibdir  = "/Users/ward/.cabal/lib/x86_64-osx-ghc-8.10.2"
datadir    = "/Users/ward/.cabal/share/x86_64-osx-ghc-8.10.2/mergeDistMats-0.1.0.0"
libexecdir = "/Users/ward/.cabal/libexec/x86_64-osx-ghc-8.10.2/mergeDistMats-0.1.0.0"
sysconfdir = "/Users/ward/.cabal/etc"

getBinDir, getLibDir, getDynLibDir, getDataDir, getLibexecDir, getSysconfDir :: IO FilePath
getBinDir = catchIO (getEnv "mergeDistMats_bindir") (\_ -> return bindir)
getLibDir = catchIO (getEnv "mergeDistMats_libdir") (\_ -> return libdir)
getDynLibDir = catchIO (getEnv "mergeDistMats_dynlibdir") (\_ -> return dynlibdir)
getDataDir = catchIO (getEnv "mergeDistMats_datadir") (\_ -> return datadir)
getLibexecDir = catchIO (getEnv "mergeDistMats_libexecdir") (\_ -> return libexecdir)
getSysconfDir = catchIO (getEnv "mergeDistMats_sysconfdir") (\_ -> return sysconfdir)

getDataFileName :: FilePath -> IO FilePath
getDataFileName name = do
  dir <- getDataDir
  return (dir ++ "/" ++ name)
