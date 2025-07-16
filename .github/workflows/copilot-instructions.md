# stimgate

This is an R package intended for eventual submission to BioConductor.
We want to identify cells that have possibly responded to stimulation, by comparing the unstimulated and stimulated tubes from the same sample.

## Code standards

### Required before each commit

- Not sure what to put here

### Development flow

## Repository structure

- `R`: Core R source code
  - 

- `src/`: Core TypeScript source code
  - `main.ts`: Main entry point and action orchestration
  - `fileValidator.ts`: Core file validation logic
  - `index.ts`: Action entrypoint that calls run()
  - `types.ts`: TypeScript type definitions
- `__tests__/`: Jest unit tests for all source files
- `dist/`: Compiled and bundled JavaScript output (generated)
- `action.yml`: GitHub Action metadata and interface definition
- `script/`: Release automation scripts
- `badges/`: Generated coverage and status badges