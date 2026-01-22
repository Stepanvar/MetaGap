# Active Tasks

## Current State
The project is currently in a mature development phase with established workflows for VCF ingestion, metadata management, and cohort analysis. The codebase appears stable with comprehensive test coverage requirements and well-defined interfaces for both the web application and command-line tools. The system supports internationalization with English and Russian interface translations, and includes production-ready deployment configurations.

## TODOs
Based on the codebase analysis, there are no explicit `@todo` or `TODO` comments in the source code. However, potential areas for enhancement based on the architecture include:
- Expanding the metadata configuration to support additional VCF header sections beyond the currently defined ones
- Improving error handling for edge cases in VCF parsing, especially for malformed files or unusual metadata structures
- Adding more sophisticated filtering capabilities for variant searches in the web interface
- Extending the VCF merging workflow to support additional file formats or more complex joint-calling scenarios
- Enhancing the documentation and examples for customizing the metadata field mappings in `metadata_fields.yaml`