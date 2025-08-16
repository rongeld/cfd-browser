"use client";
import React, { useState, useEffect, useRef } from "react";
import {
  Play,
  Pause,
  SkipForward,
  RotateCcw,
  ArrowRight,
  Droplets,
  Wind,
  Gauge,
  Zap,
  Target,
  Calculator,
} from "lucide-react";

interface SimplifiedPDEVisualizerProps {
  gridWidth?: number;
  gridHeight?: number;
}

const SimplifiedPDEVisualizer = ({
  gridWidth = 6,
  gridHeight = 4,
}: SimplifiedPDEVisualizerProps) => {
  const [currentStep, setCurrentStep] = useState(0);
  const [isAnimating, setIsAnimating] = useState(false);
  const [selectedCell, setSelectedCell] = useState<{ x: number; y: number } | null>(null);
  const [showCalculations, setShowCalculations] = useState(false);
  const [speed, setSpeed] = useState(2000); // Animation speed in ms
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [animationFrame, setAnimationFrame] = useState(0);

  // Simple 2D grid state with more detailed tracking
  const [grid, setGrid] = useState(() => {
    const size = gridWidth * gridHeight;
    return {
      velocityX: new Array(size).fill(0),
      velocityY: new Array(size).fill(0),
      pressure: new Array(size).fill(0),
      divergence: new Array(size).fill(0),
      obstacles: new Array(size).fill(false),
      // Track intermediate calculations
      flowIn: new Array(size).fill(0),
      flowOut: new Array(size).fill(0),
      pressureGradX: new Array(size).fill(0),
      pressureGradY: new Array(size).fill(0),
    };
  });

  // Initialize with simple setup
  useEffect(() => {
    const size = gridWidth * gridHeight;
    const newGrid = {
      velocityX: new Array(size).fill(0),
      velocityY: new Array(size).fill(0),
      pressure: new Array(size).fill(0),
      divergence: new Array(size).fill(0),
      obstacles: new Array(size).fill(false),
      flowIn: new Array(size).fill(0),
      flowOut: new Array(size).fill(0),
      pressureGradX: new Array(size).fill(0),
      pressureGradY: new Array(size).fill(0),
    };

    // Add a simple obstacle in the middle
    const centerX = Math.floor(gridWidth / 2);
    const centerY = Math.floor(gridHeight / 2);
    const obstacleIdx = centerY * gridWidth + centerX;
    newGrid.obstacles[obstacleIdx] = true;

    // Add initial wind from left
    for (let j = 0; j < gridHeight; j++) {
      for (let i = 0; i < 2; i++) {
        // First 2 columns
        const idx = j * gridWidth + i;
        if (!newGrid.obstacles[idx]) {
          newGrid.velocityX[idx] = 1.5; // Wind speed
          newGrid.velocityY[idx] = 0;
        }
      }
    }

    setGrid(newGrid);
  }, [gridWidth, gridHeight]);

  // Handle canvas click to select a cell
  const handleCanvasClick = (event: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const rect = canvas.getBoundingClientRect();
    const clickX = event.clientX - rect.left;
    const clickY = event.clientY - rect.top;
    
    const cellWidth = canvas.width / gridWidth;
    const cellHeight = canvas.height / gridHeight;
    
    const gridX = Math.floor(clickX / cellWidth);
    const gridY = Math.floor(clickY / cellHeight);
    
    if (gridX >= 0 && gridX < gridWidth && gridY >= 0 && gridY < gridHeight) {
      setSelectedCell({ x: gridX, y: gridY });
      setShowCalculations(true);
    }
  };

  // Calculate actual PDE values for a specific cell
  const calculateCellPDE = (x: number, y: number) => {
    const idx = y * gridWidth + x;
    
    // Divergence calculation: ∇·v = ∂vx/∂x + ∂vy/∂y
    const leftIdx = x > 0 ? idx - 1 : idx;
    const rightIdx = x < gridWidth - 1 ? idx + 1 : idx;
    const topIdx = y > 0 ? idx - gridWidth : idx;
    const bottomIdx = y < gridHeight - 1 ? idx + gridWidth : idx;
    
    const dvx_dx = (grid.velocityX[rightIdx] - grid.velocityX[leftIdx]) / 2;
    const dvy_dy = (grid.velocityY[bottomIdx] - grid.velocityY[topIdx]) / 2;
    const divergence = dvx_dx + dvy_dy;
    
    // Pressure gradients: ∇p = (∂p/∂x, ∂p/∂y)
    const dpx_dx = (grid.pressure[rightIdx] - grid.pressure[leftIdx]) / 2;
    const dpy_dy = (grid.pressure[bottomIdx] - grid.pressure[topIdx]) / 2;
    
    // Flow in/out calculation
    const flowInX = Math.max(0, grid.velocityX[leftIdx]);
    const flowOutX = Math.max(0, grid.velocityX[idx]);
    const flowInY = Math.max(0, grid.velocityY[topIdx]);
    const flowOutY = Math.max(0, grid.velocityY[idx]);
    
    const totalFlowIn = flowInX + flowInY;
    const totalFlowOut = flowOutX + flowOutY;
    
    return {
      position: { x, y },
      velocity: { x: grid.velocityX[idx], y: grid.velocityY[idx] },
      pressure: grid.pressure[idx],
      divergence,
      pressureGrad: { x: dpx_dx, y: dpy_dy },
      flow: { in: totalFlowIn, out: totalFlowOut, net: totalFlowIn - totalFlowOut },
      neighbors: {
        left: grid.velocityX[leftIdx],
        right: grid.velocityX[rightIdx],
        top: grid.velocityY[topIdx],
        bottom: grid.velocityY[bottomIdx]
      }
    };
  };

  // Teaching steps with simple explanations
  const teachingSteps = [
    {
      id: 0,
      title: "🌬️ Step 1: Wind Flows In",
      explanation:
        "Air enters from the left side. Each cell has a velocity (speed and direction).",
      highlight: "velocity",
      action: "We set the wind speed to 2 units/second flowing rightward.",
      visual: "Blue arrows show wind direction and speed.",
    },
    {
      id: 1,
      title: "🧱 Step 2: Wind Hits Obstacle",
      explanation:
        "When wind hits the obstacle, it can't flow through. What happens?",
      highlight: "obstacle",
      action:
        "The wind must flow around the obstacle - it can't go through solid objects!",
      visual: "Gray squares are solid obstacles.",
    },
    {
      id: 2,
      title: "📊 Step 3: Check Mass Conservation",
      explanation:
        "Air can't just disappear! If more air flows in than flows out, we have a problem.",
      highlight: "divergence",
      action: "We calculate: Is more air flowing IN or OUT of each cell?",
      visual: "Red = too much air coming in, Blue = too much air going out.",
    },
    {
      id: 3,
      title: "💨 Step 4: Create Pressure",
      explanation:
        "When air piles up (can't escape), pressure increases. High pressure pushes air away.",
      highlight: "pressure",
      action:
        "We solve: Where does pressure need to be high/low to balance the flow?",
      visual: "Red = high pressure, Blue = low pressure.",
    },
    {
      id: 4,
      title: "⚡ Step 5: Pressure Pushes Flow",
      explanation:
        "High pressure pushes air toward low pressure, fixing the 'too much air' problem.",
      highlight: "projection",
      action:
        "We adjust wind speed: faster away from high pressure, slower toward low pressure.",
      visual: "Arrows change to balance the flow.",
    },
    {
      id: 5,
      title: "🌊 Step 6: Transport the Air",
      explanation:
        "Air carries itself along - fast-moving air travels further in the next time step.",
      highlight: "advection",
      action:
        "We move the air: each bit of air travels according to its current speed.",
      visual: "Wind patterns move and shift position.",
    },
    {
      id: 6,
      title: "🍯 Step 7: Add Stickiness (Viscosity)",
      explanation:
        "Real air is slightly sticky - fast air slows down nearby slow air.",
      highlight: "viscosity",
      action:
        "We smooth the flow: neighboring cells influence each other's speed.",
      visual: "Sharp changes in speed become smoother.",
    },
    {
      id: 7,
      title: "🔄 Step 8: Repeat!",
      explanation:
        "This happens 60 times per second! Each cycle moves the simulation forward in time.",
      highlight: "complete",
      action: "Go back to step 1 and repeat for the next moment in time.",
      visual: "The cycle continues, creating flowing motion.",
    },
  ];

  // Simulate the teaching step
  const simulateStep = (step: number) => {
    const newGrid = { ...grid };

    switch (step) {
      case 0: // Initial wind
        for (let j = 0; j < gridHeight; j++) {
          for (let i = 0; i < 3; i++) {
            const idx = j * gridWidth + i;
            if (!newGrid.obstacles[idx]) {
              newGrid.velocityX[idx] = 2;
              newGrid.velocityY[idx] = 0;
            }
          }
        }
        break;

      case 1: // Show obstacle blocking
        // Wind starts to go around obstacle
        const centerX = Math.floor(gridWidth / 2);
        const centerY = Math.floor(gridHeight / 2);
        // Add some vertical velocity around obstacle
        if (centerY > 0) {
          const aboveIdx = (centerY - 1) * gridWidth + centerX;
          newGrid.velocityY[aboveIdx] = -0.5;
        }
        if (centerY < gridHeight - 1) {
          const belowIdx = (centerY + 1) * gridWidth + centerX;
          newGrid.velocityY[belowIdx] = 0.5;
        }
        break;

      case 2: // Calculate divergence (simplified visualization)
        // Show where too much air is coming in (red) or going out (blue)
        for (let i = 1; i < gridWidth - 1; i++) {
          for (let j = 1; j < gridHeight - 1; j++) {
            const idx = j * gridWidth + i;
            if (!newGrid.obstacles[idx]) {
              const flowIn =
                newGrid.velocityX[idx - 1] + newGrid.velocityY[idx - gridWidth];
              const flowOut = newGrid.velocityX[idx] + newGrid.velocityY[idx];
              const divergence = flowOut - flowIn;
              
              // Store divergence and flow calculations for visualization
              newGrid.divergence[idx] = divergence;
              newGrid.flowIn[idx] = flowIn;
              newGrid.flowOut[idx] = flowOut;
              
              // Store divergence in pressure field for visualization
              newGrid.pressure[idx] = divergence;
            }
          }
        }
        break;

      case 3: // Create pressure
        // High pressure where too much air is coming in
        for (let i = 0; i < gridWidth * gridHeight; i++) {
          if (!newGrid.obstacles[i]) {
            newGrid.pressure[i] = Math.random() * 0.5 - 0.25; // Random pressure for demo
          }
        }
        break;

      case 4: // Pressure affects velocity
        for (let i = 1; i < gridWidth - 1; i++) {
          for (let j = 1; j < gridHeight - 1; j++) {
            const idx = j * gridWidth + i;
            if (!newGrid.obstacles[idx]) {
              // Calculate pressure gradient
              const pressureGradX =
                newGrid.pressure[idx + 1] - newGrid.pressure[idx - 1];
              const pressureGradY =
                newGrid.pressure[idx + gridWidth] -
                newGrid.pressure[idx - gridWidth];
              
              // Store pressure gradients for visualization  
              newGrid.pressureGradX[idx] = pressureGradX;
              newGrid.pressureGradY[idx] = pressureGradY;
              
              // Apply pressure gradient to velocity
              newGrid.velocityX[idx] -= pressureGradX * 0.1;
              newGrid.velocityY[idx] -= pressureGradY * 0.1;
            }
          }
        }
        break;

      case 5: // Advection (simplified)
        // Shift velocities slightly to show transport
        const tempVelX = [...newGrid.velocityX];
        const tempVelY = [...newGrid.velocityY];
        for (let i = 1; i < gridWidth - 1; i++) {
          for (let j = 1; j < gridHeight - 1; j++) {
            const idx = j * gridWidth + i;
            if (!newGrid.obstacles[idx]) {
              // Simple advection approximation
              newGrid.velocityX[idx] =
                tempVelX[idx] * 0.98 + tempVelX[idx - 1] * 0.02;
              newGrid.velocityY[idx] =
                tempVelY[idx] * 0.98 + tempVelY[idx - gridWidth] * 0.02;
            }
          }
        }
        break;

      case 6: // Viscosity (smoothing)
        const smoothVelX = [...newGrid.velocityX];
        const smoothVelY = [...newGrid.velocityY];
        for (let i = 1; i < gridWidth - 1; i++) {
          for (let j = 1; j < gridHeight - 1; j++) {
            const idx = j * gridWidth + i;
            if (!newGrid.obstacles[idx]) {
              // Simple smoothing
              const avgX =
                (smoothVelX[idx - 1] +
                  smoothVelX[idx + 1] +
                  smoothVelX[idx - gridWidth] +
                  smoothVelX[idx + gridWidth]) /
                4;
              const avgY =
                (smoothVelY[idx - 1] +
                  smoothVelY[idx + 1] +
                  smoothVelY[idx - gridWidth] +
                  smoothVelY[idx + gridWidth]) /
                4;
              newGrid.velocityX[idx] = smoothVelX[idx] * 0.9 + avgX * 0.1;
              newGrid.velocityY[idx] = smoothVelY[idx] * 0.9 + avgY * 0.1;
            }
          }
        }
        break;

      case 7: // Complete cycle
        // Reset for next cycle
        break;
    }

    setGrid(newGrid);
  };

  // Draw the grid
  const draw = () => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const cellWidth = canvas.width / gridWidth;
    const cellHeight = canvas.height / gridHeight;

    ctx.clearRect(0, 0, canvas.width, canvas.height);

    const step = teachingSteps[currentStep];

    // Draw grid
    ctx.strokeStyle = "#ddd";
    ctx.lineWidth = 1;
    for (let i = 0; i <= gridWidth; i++) {
      ctx.beginPath();
      ctx.moveTo(i * cellWidth, 0);
      ctx.lineTo(i * cellWidth, canvas.height);
      ctx.stroke();
    }
    for (let j = 0; j <= gridHeight; j++) {
      ctx.beginPath();
      ctx.moveTo(0, j * cellHeight);
      ctx.lineTo(canvas.width, j * cellHeight);
      ctx.stroke();
    }

    // Draw based on current step
    for (let i = 0; i < gridWidth; i++) {
      for (let j = 0; j < gridHeight; j++) {
        const idx = j * gridWidth + i;
        const x = i * cellWidth;
        const y = j * cellHeight;

        // Draw obstacles
        if (grid.obstacles[idx]) {
          ctx.fillStyle = "#666";
          ctx.fillRect(x + 2, y + 2, cellWidth - 4, cellHeight - 4);
          continue;
        }

        // Draw based on what we're highlighting
        if (
          step.highlight === "velocity" ||
          step.highlight === "projection" ||
          step.highlight === "advection" ||
          step.highlight === "viscosity"
        ) {
          // Draw velocity arrows
          const vx = grid.velocityX[idx];
          const vy = grid.velocityY[idx];
          const speed = Math.sqrt(vx * vx + vy * vy);

          if (speed > 0.1) {
            const centerX = x + cellWidth / 2;
            const centerY = y + cellHeight / 2;
            const scale = Math.min(cellWidth, cellHeight) * 0.3;
            const arrowX = centerX + vx * scale;
            const arrowY = centerY + vy * scale;

            // Color by speed
            const intensity = Math.min(255, speed * 100);
            ctx.strokeStyle = `rgb(0, 100, ${intensity})`;
            ctx.fillStyle = `rgb(0, 100, ${intensity})`;
            ctx.lineWidth = 2;

            // Draw arrow
            ctx.beginPath();
            ctx.moveTo(centerX, centerY);
            ctx.lineTo(arrowX, arrowY);
            ctx.stroke();

            // Arrow head
            const angle = Math.atan2(vy, vx);
            ctx.beginPath();
            ctx.moveTo(arrowX, arrowY);
            ctx.lineTo(
              arrowX - 5 * Math.cos(angle - 0.5),
              arrowY - 5 * Math.sin(angle - 0.5)
            );
            ctx.lineTo(
              arrowX - 5 * Math.cos(angle + 0.5),
              arrowY - 5 * Math.sin(angle + 0.5)
            );
            ctx.closePath();
            ctx.fill();
          }
        } else if (
          step.highlight === "pressure" ||
          step.highlight === "divergence"
        ) {
          // Draw pressure/divergence field
          const value = grid.pressure[idx];
          const intensity = Math.min(255, Math.abs(value) * 200);

          if (value > 0) {
            ctx.fillStyle = `rgba(255, 0, 0, ${intensity / 255})`;
          } else {
            ctx.fillStyle = `rgba(0, 0, 255, ${intensity / 255})`;
          }
          ctx.fillRect(x + 2, y + 2, cellWidth - 4, cellHeight - 4);
        }
      }
    }

    // Highlight selected cell
    if (selectedCell) {
      const cellWidth = canvas.width / gridWidth;
      const cellHeight = canvas.height / gridHeight;
      const x = selectedCell.x * cellWidth;
      const y = selectedCell.y * cellHeight;
      
      ctx.strokeStyle = "#ff0000";
      ctx.lineWidth = 3;
      ctx.setLineDash([5, 5]);
      ctx.strokeRect(x, y, cellWidth, cellHeight);
      ctx.setLineDash([]);
      
      // Add a small label
      ctx.fillStyle = "#ff0000";
      ctx.font = "bold 12px Arial";
      ctx.fillText(`(${selectedCell.x},${selectedCell.y})`, x + 2, y + 15);
    }

    // Draw legend
    ctx.fillStyle = "black";
    ctx.font = "bold 14px Arial";
    ctx.fillText(step.visual, 10, canvas.height - 10);

    // Add instruction text
    if (!selectedCell) {
      ctx.fillStyle = "#666";
      ctx.font = "12px Arial";
      ctx.fillText("Click on any cell to see detailed PDE calculations", 10, 20);
    }
  };

  useEffect(() => {
    draw();
  }, [grid, currentStep, selectedCell]);

  useEffect(() => {
    if (isAnimating) {
      const interval = setInterval(() => {
        setCurrentStep((prev) => {
          const next = (prev + 1) % teachingSteps.length;
          simulateStep(next);
          return next;
        });
      }, speed);

      return () => clearInterval(interval);
    }
  }, [isAnimating, speed]);

  const handleStepForward = () => {
    const next = (currentStep + 1) % teachingSteps.length;
    setCurrentStep(next);
    simulateStep(next);
  };

  const handleReset = () => {
    setCurrentStep(0);
    setIsAnimating(false);
    // Reset grid
    useEffect(() => {
      const newGrid = {
        velocityX: new Array(gridWidth * gridHeight).fill(0),
        velocityY: new Array(gridWidth * gridHeight).fill(0),
        pressure: new Array(gridWidth * gridHeight).fill(0),
        obstacles: new Array(gridWidth * gridHeight).fill(false),
        divergence: new Array(gridWidth * gridHeight).fill(0),
        flowIn: new Array(gridWidth * gridHeight).fill(0),
        flowOut: new Array(gridWidth * gridHeight).fill(0),
        pressureGradX: new Array(gridWidth * gridHeight).fill(0),
        pressureGradY: new Array(gridWidth * gridHeight).fill(0),
      };

      // Add obstacle
      const centerX = Math.floor(gridWidth / 2);
      const centerY = Math.floor(gridHeight / 2);
      const obstacleIdx = centerY * gridWidth + centerX;
      newGrid.obstacles[obstacleIdx] = true;

      setGrid(newGrid);
    }, []);
  };

  const step = teachingSteps[currentStep];

  return (
    <div className="space-y-6">
      <div className="text-center">
        <h3 className="text-xl font-bold mb-2">
          How Does a Computer Solve Fluid Flow?
        </h3>
        <p className="text-sm text-gray-600 dark:text-gray-400">
          Let's break it down into simple steps that anyone can understand!
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
          {isAnimating ? "Pause" : "Auto Play"}
        </button>

        <button
          onClick={handleStepForward}
          className="px-4 py-2 bg-blue-500 hover:bg-blue-600 text-white rounded-lg flex items-center gap-2"
        >
          <SkipForward className="w-4 h-4" />
          Next Step
        </button>

        <button
          onClick={handleReset}
          className="px-4 py-2 bg-gray-500 hover:bg-gray-600 text-white rounded-lg flex items-center gap-2"
        >
          <RotateCcw className="w-4 h-4" />
          Reset
        </button>
      </div>

      {/* Speed Control */}
      <div className="flex justify-center items-center gap-4">
        <span className="text-sm">Speed:</span>
        <input
          type="range"
          min="500"
          max="3000"
          step="100"
          value={speed}
          onChange={(e) => setSpeed(parseInt(e.target.value))}
          className="w-32"
        />
        <span className="text-sm">{(speed / 1000).toFixed(1)}s</span>
      </div>

      {/* Step Progress */}
      <div className="flex justify-center">
        <div className="flex gap-2">
          {teachingSteps.map((_, index) => (
            <div
              key={index}
              className={`w-3 h-3 rounded-full ${
                index === currentStep ? "bg-blue-500" : "bg-gray-300"
              }`}
            />
          ))}
        </div>
      </div>

      {/* Visualization */}
      <div className="flex justify-center">
        <canvas
          ref={canvasRef}
          width={480}
          height={360}
          onClick={handleCanvasClick}
          className="border-2 border-gray-300 dark:border-gray-600 rounded-lg bg-white cursor-crosshair"
        />
      </div>

      {/* Current Step Explanation */}
      <div className="bg-blue-50 dark:bg-blue-900/20 p-6 rounded-lg max-w-2xl mx-auto">
        <div className="flex items-center gap-3 mb-3">
          <div className="text-2xl">{step.title.split(" ")[0]}</div>
          <h4 className="text-lg font-semibold">{step.title.slice(2)}</h4>
        </div>

        <div className="space-y-3">
          <div className="bg-white dark:bg-gray-800 p-3 rounded">
            <strong>What's happening:</strong> {step.explanation}
          </div>

          <div className="bg-green-50 dark:bg-green-900/20 p-3 rounded">
            <strong>Computer does:</strong> {step.action}
          </div>

          <div className="bg-yellow-50 dark:bg-yellow-900/20 p-3 rounded">
            <strong>You see:</strong> {step.visual}
          </div>
        </div>
      </div>

      {/* Simple Summary */}
      <div className="bg-gray-50 dark:bg-gray-800 p-4 rounded-lg max-w-2xl mx-auto">
        <h4 className="font-semibold mb-2">The Big Picture:</h4>
        <p className="text-sm">
          The computer breaks time into tiny slices (1/60th second). For each
          slice, it:
          <br />
          1. Figures out where air wants to go
          <br />
          2. Makes sure air doesn't pile up or disappear
          <br />
          3. Lets air flow naturally
          <br />
          4. Adds realistic effects like stickiness
          <br />
          <br />
          Do this 60 times per second = smooth flowing animation! 🎬
        </p>
      </div>

      {/* Detailed Cell Calculations */}
      {selectedCell && showCalculations && (
        <div className="bg-green-50 dark:bg-green-900/20 p-6 rounded-lg max-w-4xl mx-auto mt-4">
          <div className="flex justify-between items-start mb-4">
            <div className="flex items-center gap-3">
              <div className="text-2xl">🔬</div>
              <h4 className="text-lg font-semibold">
                Cell ({selectedCell.x}, {selectedCell.y}) - Detailed PDE Calculations
              </h4>
            </div>
            <button 
              onClick={() => setShowCalculations(false)}
              className="text-gray-500 hover:text-gray-700 text-xl"
            >
              ×
            </button>
          </div>

          {(() => {
            const calc = calculateCellPDE(selectedCell.x, selectedCell.y);
            return (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                {/* Current State */}
                <div className="space-y-3">
                  <h5 className="font-semibold text-blue-700 dark:text-blue-300">📍 Current State</h5>
                  <div className="bg-white dark:bg-gray-800 p-3 rounded">
                    <div className="text-sm space-y-1">
                      <div><strong>Position:</strong> Grid cell ({calc.position.x}, {calc.position.y})</div>
                      <div><strong>Velocity X:</strong> {calc.velocity.x.toFixed(3)} m/s</div>
                      <div><strong>Velocity Y:</strong> {calc.velocity.y.toFixed(3)} m/s</div>
                      <div><strong>Speed:</strong> {Math.sqrt(calc.velocity.x**2 + calc.velocity.y**2).toFixed(3)} m/s</div>
                      <div><strong>Pressure:</strong> {calc.pressure.toFixed(3)} Pa</div>
                    </div>
                  </div>
                </div>

                {/* PDE Calculations */}
                <div className="space-y-3">
                  <h5 className="font-semibold text-purple-700 dark:text-purple-300">🧮 PDE Math</h5>
                  <div className="bg-white dark:bg-gray-800 p-3 rounded">
                    <div className="text-sm space-y-2">
                      <div>
                        <strong>Divergence (∇·v):</strong> {calc.divergence.toFixed(4)}
                        <div className="text-xs text-gray-600 ml-4">
                          = ∂vₓ/∂x + ∂vᵧ/∂y<br/>
                          = {(calc.neighbors.right - calc.neighbors.left)/2} + {(calc.neighbors.bottom - calc.neighbors.top)/2}
                        </div>
                      </div>
                      <div>
                        <strong>Pressure Gradient:</strong> 
                        <div className="text-xs text-gray-600 ml-4">
                          ∇p = ({calc.pressureGrad.x.toFixed(4)}, {calc.pressureGrad.y.toFixed(4)})
                        </div>
                      </div>
                    </div>
                  </div>
                </div>

                {/* Flow Analysis */}
                <div className="space-y-3">
                  <h5 className="font-semibold text-green-700 dark:text-green-300">🌊 Flow Analysis</h5>
                  <div className="bg-white dark:bg-gray-800 p-3 rounded">
                    <div className="text-sm space-y-1">
                      <div><strong>Flow In:</strong> {calc.flow.in.toFixed(4)} m³/s</div>
                      <div><strong>Flow Out:</strong> {calc.flow.out.toFixed(4)} m³/s</div>
                      <div><strong>Net Flow:</strong> {calc.flow.net.toFixed(4)} m³/s</div>
                      <div className="text-xs text-gray-600 mt-2">
                        {Math.abs(calc.flow.net) < 0.001 ? 
                          "✅ Balanced flow (conservation satisfied)" : 
                          "⚠️ Imbalanced flow (pressure correction needed)"
                        }
                      </div>
                    </div>
                  </div>
                </div>

                {/* Neighbor Values */}
                <div className="space-y-3">
                  <h5 className="font-semibold text-orange-700 dark:text-orange-300">🔄 Neighbor Interaction</h5>
                  <div className="bg-white dark:bg-gray-800 p-3 rounded">
                    <div className="text-sm">
                      <div className="grid grid-cols-3 gap-1 text-center mb-2">
                        <div></div>
                        <div className="bg-blue-100 p-1 rounded">↑ {calc.neighbors.top.toFixed(2)}</div>
                        <div></div>
                        <div className="bg-blue-100 p-1 rounded">← {calc.neighbors.left.toFixed(2)}</div>
                        <div className="bg-yellow-200 p-1 rounded font-bold">
                          ({selectedCell.x},{selectedCell.y})
                        </div>
                        <div className="bg-blue-100 p-1 rounded">→ {calc.neighbors.right.toFixed(2)}</div>
                        <div></div>
                        <div className="bg-blue-100 p-1 rounded">↓ {calc.neighbors.bottom.toFixed(2)}</div>
                        <div></div>
                      </div>
                      <div className="text-xs text-gray-600 text-center">
                        Velocities of neighboring cells
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            );
          })()}

          <div className="mt-4 p-3 bg-blue-100 dark:bg-blue-900/30 rounded text-sm">
            <strong>💡 What this means:</strong> These calculations happen simultaneously for all {gridWidth * gridHeight} cells every simulation step. 
            The computer uses these values to update velocities and pressures according to the Navier-Stokes equations!
          </div>
        </div>
      )}
    </div>
  );
};

export default SimplifiedPDEVisualizer;
