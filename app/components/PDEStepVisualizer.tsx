"use client";
import React, { useState, useEffect, useRef } from "react";
import { Play, Pause, Eye, Calculator, TrendingUp } from "lucide-react";

interface PDEVisualizerProps {
  gridWidth: number;
  gridHeight: number;
  simulationData?: any;
  onStepUpdate?: (step: string, data: any) => void;
}

const PDEStepVisualizer = ({
  gridWidth,
  gridHeight,
  simulationData,
  onStepUpdate,
}: PDEVisualizerProps) => {
  const [activeStep, setActiveStep] = useState<
    "divergence" | "pressure" | "projection" | "advection" | "viscosity"
  >("divergence");
  const [iterationCount, setIterationCount] = useState(0);
  const [isAnimating, setIsAnimating] = useState(false);
  const [showValues, setShowValues] = useState(true);
  const canvasRef = useRef<HTMLCanvasElement>(null);

  // Mock data for demonstration - in real implementation this would come from the actual simulation
  const generateMockData = (step: string, iteration: number = 0) => {
    const size = Math.min(gridWidth, gridHeight);
    const data = new Array(size * size);

    for (let i = 0; i < size; i++) {
      for (let j = 0; j < size; j++) {
        const idx = i * size + j;
        switch (step) {
          case "divergence":
            data[idx] = Math.sin(i * 0.3) * Math.cos(j * 0.3) * 0.5;
            break;
          case "pressure":
            // Simulate Jacobi iterations converging
            const convergence = Math.max(0, 1 - iteration * 0.1);
            data[idx] =
              Math.sin(i * 0.2 + iteration * 0.1) *
              Math.cos(j * 0.2) *
              convergence;
            break;
          case "projection":
            data[idx] = Math.cos(i * 0.25) * Math.sin(j * 0.25) * 0.3;
            break;
          case "advection":
            data[idx] =
              Math.sin((i + iteration * 0.5) * 0.2) * Math.cos(j * 0.2) * 0.4;
            break;
          case "viscosity":
            data[idx] =
              Math.exp(
                -((i - size / 2) ** 2 + (j - size / 2) ** 2) / (size * 0.3)
              ) * Math.sin(iteration * 0.3);
            break;
          default:
            data[idx] = 0;
        }
      }
    }
    return data;
  };

  const drawField = (
    canvas: HTMLCanvasElement,
    data: number[],
    title: string
  ) => {
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const size = Math.sqrt(data.length);
    const cellSize = Math.min(canvas.width / size, canvas.height / size);

    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Find min/max for color scaling
    const min = Math.min(...data);
    const max = Math.max(...data);
    const range = max - min || 1;

    for (let i = 0; i < size; i++) {
      for (let j = 0; j < size; j++) {
        const idx = i * size + j;
        const value = data[idx];
        const normalizedValue = (value - min) / range;

        // Color mapping: blue for negative, red for positive
        const r = value > 0 ? Math.floor(normalizedValue * 255) : 0;
        const b = value < 0 ? Math.floor((1 - normalizedValue) * 255) : 0;
        const g = Math.floor((1 - Math.abs(normalizedValue)) * 100);

        ctx.fillStyle = `rgb(${r}, ${g}, ${b})`;
        ctx.fillRect(j * cellSize, i * cellSize, cellSize - 1, cellSize - 1);

        // Show values if enabled and cell is large enough
        if (showValues && cellSize > 20) {
          ctx.fillStyle = "white";
          ctx.font = `${cellSize * 0.3}px monospace`;
          ctx.textAlign = "center";
          ctx.fillText(
            value.toFixed(2),
            j * cellSize + cellSize / 2,
            i * cellSize + cellSize / 2 + cellSize * 0.1
          );
        }
      }
    }

    // Title
    ctx.fillStyle = "black";
    ctx.font = "bold 16px Arial";
    ctx.textAlign = "left";
    ctx.fillText(title, 10, 20);
  };

  useEffect(() => {
    if (!isAnimating) return;

    const interval = setInterval(() => {
      setIterationCount((prev) => {
        const next = (prev + 1) % 20;

        // Cycle through steps
        if (next === 0) {
          const steps: (typeof activeStep)[] = [
            "divergence",
            "pressure",
            "projection",
            "advection",
            "viscosity",
          ];
          const currentIndex = steps.indexOf(activeStep);
          const nextStep = steps[(currentIndex + 1) % steps.length];
          setActiveStep(nextStep);
        }

        return next;
      });
    }, 200);

    return () => clearInterval(interval);
  }, [isAnimating, activeStep]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const data = generateMockData(activeStep, iterationCount);
    drawField(
      canvas,
      data,
      `${activeStep.toUpperCase()} - Iteration ${iterationCount}`
    );
  }, [activeStep, iterationCount, showValues, gridWidth, gridHeight]);

  const stepDescriptions = {
    divergence: {
      title: "1. Divergence Calculation",
      formula: "∇·u = ∂u/∂x + ∂v/∂y",
      description:
        "Calculate how much the velocity field is 'spreading out' at each point. Non-zero divergence violates mass conservation.",
      color: "bg-red-50 dark:bg-red-900/20",
    },
    pressure: {
      title: "2. Pressure Poisson Equation",
      formula: "∇²p = -ρ∇·u",
      description:
        "Solve for pressure using Jacobi iterations. This is the most computationally expensive step.",
      color: "bg-blue-50 dark:bg-blue-900/20",
    },
    projection: {
      title: "3. Velocity Projection",
      formula: "u_new = u - ∇p/ρ",
      description:
        "Subtract pressure gradient from velocity to make the flow incompressible (divergence-free).",
      color: "bg-green-50 dark:bg-green-900/20",
    },
    advection: {
      title: "4. Semi-Lagrangian Advection",
      formula: "∂u/∂t + (u·∇)u = 0",
      description:
        "Transport velocity by following particles backwards in time. This handles the nonlinear convection term.",
      color: "bg-purple-50 dark:bg-purple-900/20",
    },
    viscosity: {
      title: "5. Viscosity Diffusion",
      formula: "∂u/∂t = ν∇²u",
      description:
        "Apply viscous diffusion using explicit finite differences. Makes the flow smoother.",
      color: "bg-yellow-50 dark:bg-yellow-900/20",
    },
  };

  return (
    <div className="space-y-6">
      <div className="text-center">
        <h3 className="text-xl font-bold mb-2">
          PDE Solution Process Visualizer
        </h3>
        <p className="text-sm text-gray-600 dark:text-gray-400">
          See how the simplified Navier-Stokes equations are solved step by step
        </p>
      </div>

      {/* Controls */}
      <div className="flex justify-center gap-4 mb-6">
        <button
          onClick={() => setIsAnimating(!isAnimating)}
          className={`px-4 py-2 rounded-lg flex items-center gap-2 ${
            isAnimating
              ? "bg-red-500 hover:bg-red-600"
              : "bg-green-500 hover:bg-green-600"
          } text-white transition-colors`}
        >
          {isAnimating ? (
            <Pause className="w-4 h-4" />
          ) : (
            <Play className="w-4 h-4" />
          )}
          {isAnimating ? "Pause" : "Play"}
        </button>

        <button
          onClick={() => setShowValues(!showValues)}
          className={`px-4 py-2 rounded-lg flex items-center gap-2 ${
            showValues
              ? "bg-blue-500 hover:bg-blue-600"
              : "bg-gray-500 hover:bg-gray-600"
          } text-white transition-colors`}
        >
          <Eye className="w-4 h-4" />
          {showValues ? "Hide Values" : "Show Values"}
        </button>
      </div>

      {/* Step Selection */}
      <div className="grid grid-cols-5 gap-2 mb-6">
        {Object.entries(stepDescriptions).map(([key, step]) => (
          <button
            key={key}
            onClick={() => setActiveStep(key as typeof activeStep)}
            className={`p-3 rounded-lg text-center transition-colors ${
              activeStep === key
                ? "bg-blue-500 text-white"
                : "bg-gray-200 dark:bg-gray-700 hover:bg-gray-300 dark:hover:bg-gray-600"
            }`}
          >
            <div className="text-xs font-semibold">
              {step.title.split(".")[1]}
            </div>
          </button>
        ))}
      </div>

      {/* Visualization Canvas */}
      <div className="flex justify-center">
        <canvas
          ref={canvasRef}
          width={400}
          height={400}
          className="border border-gray-300 dark:border-gray-600 rounded-lg"
        />
      </div>

      {/* Current Step Description */}
      <div className={`p-4 rounded-lg ${stepDescriptions[activeStep].color}`}>
        <div className="flex items-center gap-2 mb-2">
          <Calculator className="w-5 h-5" />
          <h4 className="font-semibold">
            {stepDescriptions[activeStep].title}
          </h4>
        </div>

        <div className="bg-white dark:bg-gray-800 p-2 rounded mb-2 font-mono text-center">
          {stepDescriptions[activeStep].formula}
        </div>

        <p className="text-sm">{stepDescriptions[activeStep].description}</p>

        {activeStep === "pressure" && (
          <div className="mt-2 text-xs text-gray-600 dark:text-gray-400">
            <strong>Jacobi Iteration {iterationCount}/20:</strong> Each
            iteration improves the pressure solution. More iterations = more
            accurate but slower.
          </div>
        )}
      </div>

      {/* Algorithm Explanation */}
      <div className="bg-gray-50 dark:bg-gray-800 p-4 rounded-lg">
        <h4 className="font-semibold mb-2 flex items-center gap-2">
          <TrendingUp className="w-4 h-4" />
          What's Actually Being Solved?
        </h4>
        <div className="text-sm space-y-2">
          <p>
            <strong>This is NOT full Navier-Stokes!</strong> It's a simplified
            incompressible flow solver using:
          </p>
          <ul className="list-disc ml-4 space-y-1">
            <li>
              <strong>Pressure Projection Method:</strong> Splits each timestep
              into advection → diffusion → projection
            </li>
            <li>
              <strong>Semi-Lagrangian Advection:</strong> Handles convection
              term (u·∇)u by particle tracing
            </li>
            <li>
              <strong>Explicit Viscosity:</strong> Simple diffusion instead of
              implicit solving
            </li>
            <li>
              <strong>Missing Terms:</strong> No body forces, simplified
              boundary conditions
            </li>
          </ul>
          <p className="mt-2 text-xs italic">
            Real CFD codes use implicit methods, higher-order discretizations,
            and sophisticated turbulence models.
          </p>
        </div>
      </div>
    </div>
  );
};

export default PDEStepVisualizer;
