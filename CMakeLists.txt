SET (src_dir ${CMAKE_CURRENT_SOURCE_DIR})

SET (FlowQA_sources
  ${src_dir}/FlowQA.cc
  ${src_dir}/FlowQA.h
  PARENT_SCOPE
) 

SET (FlowQA_data
  ${src_dir}/FlowQA.xml
  ${src_dir}/FlowQA.xsd
  PARENT_SCOPE
)
