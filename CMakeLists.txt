
cmake_minimum_required(VERSION 3.9)
cmake_policy(VERSION 3.9)
set(CMAKE_CXX_STANDARD 14)

project(alterpcb-tlinesim VERSION 0.0.0)

include(FeatureSummary)

add_subdirectory(src)

feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)
