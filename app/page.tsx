"use client";
import React, { useRef, useEffect, useState, useCallback } from "react";
import {
  Play,
  Pause,
  RotateCcw,
  Settings,
  Wind,
  Droplets,
  Brush,
  Eraser,
} from "lucide-react";

const CFDSimulator = () => {
  const canvasRef = useRef(null);
  const animationRef = useRef(null);
  const [isRunning, setIsRunning] = useState(false);
  const [drawingMode, setDrawingMode] = useState("none"); // 'draw', 'erase', 'bezier', 'none'
  const [isDrawing, setIsDrawing] = useState(false);
  const [brushSize, setBrushSize] = useState(10);
  const [settings, setSettings] = useState({
    viscosity: 0.02,
    density: 1.0,
    windSpeed: 5.0,
    particleCount: 1000,
    showVelocityField: true,
    showPressureField: false,
    showParticles: true,
    obstacleSize: 50,
    arrowSize: 1.0,
  });

  // Use ref to store current settings for simulation functions
  const settingsRef = useRef(settings);
  useEffect(() => {
    settingsRef.current = settings;
  }, [settings]);

  const [visualizationMode, setVisualizationMode] = useState("standard"); // 'standard' or 'pressure'
  const [quality, setQuality] = useState("medium"); // 'low', 'medium', 'high', 'ultra'
  const [bezierPoints, setBezierPoints] = useState([]);
  const [currentBezier, setCurrentBezier] = useState([]);

  // Grid dimensions based on quality
  const getGridDimensions = () => {
    const qualitySettings = {
      low: { width: 80, height: 40, cellSize: 10 },
      medium: { width: 120, height: 60, cellSize: 8 },
      high: { width: 160, height: 80, cellSize: 6 },
      ultra: { width: 200, height: 100, cellSize: 5 },
    };
    return qualitySettings[quality];
  };

  const {
    width: GRID_WIDTH,
    height: GRID_HEIGHT,
    cellSize: CELL_SIZE,
  } = getGridDimensions();

  // Simulation state
  const simulationState = useRef({
    velocityX: new Array(GRID_WIDTH * GRID_HEIGHT).fill(0),
    velocityY: new Array(GRID_WIDTH * GRID_HEIGHT).fill(0),
    pressure: new Array(GRID_WIDTH * GRID_HEIGHT).fill(0),
    obstacles: new Array(GRID_WIDTH * GRID_HEIGHT).fill(false),
    particles: [],
    currentGridWidth: GRID_WIDTH,
    currentGridHeight: GRID_HEIGHT,
  });

  // Reinitialize simulation when quality changes
  const reinitializeSimulation = useCallback(() => {
    const gridSize = GRID_WIDTH * GRID_HEIGHT;
    simulationState.current = {
      velocityX: new Array(gridSize).fill(0),
      velocityY: new Array(gridSize).fill(0),
      pressure: new Array(gridSize).fill(0),
      obstacles: new Array(gridSize).fill(false),
      particles: [],
      currentGridWidth: GRID_WIDTH,
      currentGridHeight: GRID_HEIGHT,
    };
  }, [GRID_WIDTH, GRID_HEIGHT]);

  // Initialize particles
  const initializeParticles = useCallback(() => {
    const particles = [];
    for (let i = 0; i < settingsRef.current.particleCount; i++) {
      particles.push({
        x: Math.random() * GRID_WIDTH * CELL_SIZE,
        y: Math.random() * GRID_HEIGHT * CELL_SIZE,
        vx: 0,
        vy: 0,
        life: 1.0,
      });
    }
    simulationState.current.particles = particles;
  }, []); // Remove dependency - will be called manually when needed

  // Initialize obstacles (start with no obstacles)
  const initializeObstacles = useCallback(() => {
    const obstacles = simulationState.current.obstacles;
    obstacles.fill(false);

    // Optional: Add default circular obstacle
    if (settingsRef.current.obstacleSize > 0) {
      const centerX = Math.floor(GRID_WIDTH * 0.3);
      const centerY = Math.floor(GRID_HEIGHT * 0.5);
      const radius = settingsRef.current.obstacleSize / CELL_SIZE;

      for (let i = 0; i < GRID_WIDTH; i++) {
        for (let j = 0; j < GRID_HEIGHT; j++) {
          const dx = i - centerX;
          const dy = j - centerY;
          if (dx * dx + dy * dy < radius * radius) {
            obstacles[j * GRID_WIDTH + i] = true;
          }
        }
      }
    }
  }, []); // Remove dependency - will be called manually when needed

  // Get grid index from coordinates
  const getIndex = (x, y) => {
    const i = Math.max(0, Math.min(GRID_WIDTH - 1, Math.floor(x)));
    const j = Math.max(0, Math.min(GRID_HEIGHT - 1, Math.floor(y)));
    return j * GRID_WIDTH + i;
  };

  // Bilinear interpolation for velocity
  const interpolateVelocity = (x, y, field) => {
    const x0 = Math.floor(x);
    const y0 = Math.floor(y);
    const x1 = Math.min(x0 + 1, GRID_WIDTH - 1);
    const y1 = Math.min(y0 + 1, GRID_HEIGHT - 1);

    const fx = x - x0;
    const fy = y - y0;

    const v00 = field[y0 * GRID_WIDTH + x0] || 0;
    const v10 = field[y0 * GRID_WIDTH + x1] || 0;
    const v01 = field[y1 * GRID_WIDTH + x0] || 0;
    const v11 = field[y1 * GRID_WIDTH + x1] || 0;

    return (
      v00 * (1 - fx) * (1 - fy) +
      v10 * fx * (1 - fy) +
      v01 * (1 - fx) * fy +
      v11 * fx * fy
    );
  };

  // Apply boundary conditions
  const applyBoundaryConditions = () => {
    const { velocityX, velocityY, obstacles } = simulationState.current;

    // Set boundary velocities
    for (let i = 0; i < GRID_WIDTH; i++) {
      for (let j = 0; j < GRID_HEIGHT; j++) {
        const idx = j * GRID_WIDTH + i;

        // Left boundary - inlet
        if (i === 0) {
          velocityX[idx] = settingsRef.current.windSpeed;
          velocityY[idx] = 0;
        }

        // Right boundary - outlet (free flow)
        if (i === GRID_WIDTH - 1) {
          velocityX[idx] = velocityX[idx - 1];
          velocityY[idx] = velocityY[idx - 1];
        }

        // Top and bottom boundaries - no slip
        if (j === 0 || j === GRID_HEIGHT - 1) {
          velocityX[idx] = 0;
          velocityY[idx] = 0;
        }

        // Obstacles - no slip
        if (obstacles[idx]) {
          velocityX[idx] = 0;
          velocityY[idx] = 0;
        }
      }
    }
  };

  // Pressure projection step
  const projectVelocity = () => {
    const { velocityX, velocityY, pressure, obstacles } =
      simulationState.current;
    const divergence = new Array(GRID_WIDTH * GRID_HEIGHT).fill(0);

    // Compute divergence
    for (let i = 1; i < GRID_WIDTH - 1; i++) {
      for (let j = 1; j < GRID_HEIGHT - 1; j++) {
        const idx = j * GRID_WIDTH + i;
        if (!obstacles[idx]) {
          const dvx = velocityX[idx + 1] - velocityX[idx - 1];
          const dvy = velocityY[idx + GRID_WIDTH] - velocityY[idx - GRID_WIDTH];
          divergence[idx] = 0.5 * (dvx + dvy);
        }
      }
    }

    // Solve pressure using Jacobi iteration
    pressure.fill(0);
    for (let iter = 0; iter < 20; iter++) {
      const newPressure = [...pressure];
      for (let i = 1; i < GRID_WIDTH - 1; i++) {
        for (let j = 1; j < GRID_HEIGHT - 1; j++) {
          const idx = j * GRID_WIDTH + i;
          if (!obstacles[idx]) {
            const neighbors =
              pressure[idx - 1] +
              pressure[idx + 1] +
              pressure[idx - GRID_WIDTH] +
              pressure[idx + GRID_WIDTH];
            newPressure[idx] = 0.25 * (neighbors - divergence[idx]);
          }
        }
      }
      pressure.splice(0, pressure.length, ...newPressure);
    }

    // Subtract pressure gradient
    for (let i = 1; i < GRID_WIDTH - 1; i++) {
      for (let j = 1; j < GRID_HEIGHT - 1; j++) {
        const idx = j * GRID_WIDTH + i;
        if (!obstacles[idx]) {
          const dpx = 0.5 * (pressure[idx + 1] - pressure[idx - 1]);
          const dpy =
            0.5 * (pressure[idx + GRID_WIDTH] - pressure[idx - GRID_WIDTH]);
          velocityX[idx] -= dpx;
          velocityY[idx] -= dpy;
        }
      }
    }
  };

  // Advection step
  const advectVelocity = () => {
    const { velocityX, velocityY } = simulationState.current;
    const newVelX = [...velocityX];
    const newVelY = [...velocityY];
    const dt = 0.016; // ~60fps

    for (let i = 1; i < GRID_WIDTH - 1; i++) {
      for (let j = 1; j < GRID_HEIGHT - 1; j++) {
        const idx = j * GRID_WIDTH + i;

        // Trace particle backwards
        const vx = velocityX[idx];
        const vy = velocityY[idx];
        const x = i - dt * vx;
        const y = j - dt * vy;

        newVelX[idx] = interpolateVelocity(x, y, velocityX);
        newVelY[idx] = interpolateVelocity(x, y, velocityY);
      }
    }

    simulationState.current.velocityX = newVelX;
    simulationState.current.velocityY = newVelY;
  };

  // Add viscosity
  const applyViscosity = () => {
    const { velocityX, velocityY, obstacles } = simulationState.current;
    const alpha = settingsRef.current.viscosity;

    for (let field of [velocityX, velocityY]) {
      const newField = [...field];
      for (let i = 1; i < GRID_WIDTH - 1; i++) {
        for (let j = 1; j < GRID_HEIGHT - 1; j++) {
          const idx = j * GRID_WIDTH + i;
          if (!obstacles[idx]) {
            const neighbors =
              field[idx - 1] +
              field[idx + 1] +
              field[idx - GRID_WIDTH] +
              field[idx + GRID_WIDTH];
            newField[idx] = field[idx] + alpha * (neighbors - 4 * field[idx]);
          }
        }
      }
      field.splice(0, field.length, ...newField);
    }
  };

  // Render function
  const render = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext("2d");
    const width = canvas.width;
    const height = canvas.height;

    // Clear canvas
    ctx.fillStyle = visualizationMode === "pressure" ? "#000011" : "#0a0a0a";
    ctx.fillRect(0, 0, width, height);

    const { velocityX, velocityY, pressure, obstacles, particles } =
      simulationState.current;

    if (visualizationMode === "pressure") {
      // Enhanced pressure visualization mode

      // Find pressure range for better scaling
      let minPressure = Infinity;
      let maxPressure = -Infinity;
      for (let i = 0; i < GRID_WIDTH * GRID_HEIGHT; i++) {
        if (!obstacles[i]) {
          minPressure = Math.min(minPressure, pressure[i]);
          maxPressure = Math.max(maxPressure, pressure[i]);
        }
      }

      const pressureRange = Math.max(
        Math.abs(minPressure),
        Math.abs(maxPressure)
      );

      // Draw enhanced pressure field
      for (let i = 0; i < GRID_WIDTH; i++) {
        for (let j = 0; j < GRID_HEIGHT; j++) {
          const idx = j * GRID_WIDTH + i;
          if (!obstacles[idx]) {
            const p = pressure[idx];
            const normalizedPressure =
              pressureRange > 0 ? p / pressureRange : 0;

            let r, g, b, alpha;
            if (normalizedPressure > 0) {
              // High pressure - red to yellow
              const intensity = Math.min(normalizedPressure * 2, 1);
              r = 255;
              g = Math.floor(intensity * 255);
              b = 0;
              alpha = 0.3 + intensity * 0.7;
            } else {
              // Low pressure - blue to cyan
              const intensity = Math.min(Math.abs(normalizedPressure) * 2, 1);
              r = 0;
              g = Math.floor(intensity * 255);
              b = 255;
              alpha = 0.3 + intensity * 0.7;
            }

            ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${alpha})`;
            ctx.fillRect(i * CELL_SIZE, j * CELL_SIZE, CELL_SIZE, CELL_SIZE);
          }
        }
      }

      // Draw pressure contour lines
      ctx.strokeStyle = "rgba(255, 255, 255, 0.3)";
      ctx.lineWidth = 1;
      const contourLevels = [-0.5, -0.25, 0, 0.25, 0.5];

      for (const level of contourLevels) {
        ctx.beginPath();
        for (let i = 1; i < GRID_WIDTH - 1; i++) {
          for (let j = 1; j < GRID_HEIGHT - 1; j++) {
            const idx = j * GRID_WIDTH + i;
            if (!obstacles[idx]) {
              const p = pressureRange > 0 ? pressure[idx] / pressureRange : 0;

              // Simple contour detection
              const neighbors = [
                pressureRange > 0 ? pressure[idx - 1] / pressureRange : 0,
                pressureRange > 0 ? pressure[idx + 1] / pressureRange : 0,
                pressureRange > 0
                  ? pressure[idx - GRID_WIDTH] / pressureRange
                  : 0,
                pressureRange > 0
                  ? pressure[idx + GRID_WIDTH] / pressureRange
                  : 0,
              ];

              const crossesLevel = neighbors.some(
                (n) => (p <= level && n >= level) || (p >= level && n <= level)
              );

              if (crossesLevel) {
                const x = i * CELL_SIZE + CELL_SIZE / 2;
                const y = j * CELL_SIZE + CELL_SIZE / 2;
                ctx.moveTo(x - 2, y);
                ctx.lineTo(x + 2, y);
                ctx.moveTo(x, y - 2);
                ctx.lineTo(x, y + 2);
              }
            }
          }
        }
        ctx.stroke();
      }

      // Minimal velocity vectors for pressure mode
      if (settingsRef.current.showVelocityField) {
        for (let i = 0; i < GRID_WIDTH; i += 4) {
          for (let j = 0; j < GRID_HEIGHT; j += 4) {
            const idx = j * GRID_WIDTH + i;
            if (!obstacles[idx]) {
              const vx = velocityX[idx];
              const vy = velocityY[idx];
              const speed = Math.sqrt(vx * vx + vy * vy);

              if (speed > 0.2) {
                const x = i * CELL_SIZE;
                const y = j * CELL_SIZE;
                const length = Math.min(
                  speed * 3 * settingsRef.current.arrowSize,
                  CELL_SIZE * 0.8
                );

                ctx.strokeStyle = `rgba(255, 255, 255, 0.6)`;
                ctx.lineWidth = 1;
                ctx.beginPath();
                ctx.moveTo(x, y);
                ctx.lineTo(x + vx * length, y + vy * length);
                ctx.stroke();

                // Arrow head
                const angle = Math.atan2(vy, vx);
                const headLength = 3 * settingsRef.current.arrowSize;
                ctx.beginPath();
                ctx.moveTo(x + vx * length, y + vy * length);
                ctx.lineTo(
                  x + vx * length - headLength * Math.cos(angle - Math.PI / 6),
                  y + vy * length - headLength * Math.sin(angle - Math.PI / 6)
                );
                ctx.moveTo(x + vx * length, y + vy * length);
                ctx.lineTo(
                  x + vx * length - headLength * Math.cos(angle + Math.PI / 6),
                  y + vy * length - headLength * Math.sin(angle + Math.PI / 6)
                );
                ctx.stroke();
              }
            }
          }
        }
      }
    } else {
      // Standard visualization mode

      // Draw velocity field
      if (settingsRef.current.showVelocityField) {
        for (let i = 0; i < GRID_WIDTH; i += 2) {
          for (let j = 0; j < GRID_HEIGHT; j += 2) {
            const idx = j * GRID_WIDTH + i;
            if (!obstacles[idx]) {
              const vx = velocityX[idx];
              const vy = velocityY[idx];
              const speed = Math.sqrt(vx * vx + vy * vy);

              if (speed > 0.1) {
                const x = i * CELL_SIZE;
                const y = j * CELL_SIZE;
                const length = Math.min(
                  speed * 5 * settingsRef.current.arrowSize,
                  CELL_SIZE
                );

                ctx.strokeStyle = `hsl(${200 + speed * 20}, 70%, 50%)`;
                ctx.lineWidth = 1;
                ctx.beginPath();
                ctx.moveTo(x, y);
                ctx.lineTo(x + vx * length, y + vy * length);
                ctx.stroke();

                // Arrow head for standard mode
                if (length > 5) {
                  const angle = Math.atan2(vy, vx);
                  const headLength = 2 * settingsRef.current.arrowSize;
                  ctx.beginPath();
                  ctx.moveTo(x + vx * length, y + vy * length);
                  ctx.lineTo(
                    x +
                      vx * length -
                      headLength * Math.cos(angle - Math.PI / 6),
                    y + vy * length - headLength * Math.sin(angle - Math.PI / 6)
                  );
                  ctx.moveTo(x + vx * length, y + vy * length);
                  ctx.lineTo(
                    x +
                      vx * length -
                      headLength * Math.cos(angle + Math.PI / 6),
                    y + vy * length - headLength * Math.sin(angle + Math.PI / 6)
                  );
                  ctx.stroke();
                }
              }
            }
          }
        }
      }

      // Draw pressure field
      if (settingsRef.current.showPressureField) {
        for (let i = 0; i < GRID_WIDTH; i++) {
          for (let j = 0; j < GRID_HEIGHT; j++) {
            const idx = j * GRID_WIDTH + i;
            if (!obstacles[idx]) {
              const p = pressure[idx];
              const intensity = Math.abs(p) * 100;
              const color =
                p > 0
                  ? `rgba(255, 0, 0, ${Math.min(intensity, 0.5)})`
                  : `rgba(0, 0, 255, ${Math.min(intensity, 0.5)})`;

              ctx.fillStyle = color;
              ctx.fillRect(i * CELL_SIZE, j * CELL_SIZE, CELL_SIZE, CELL_SIZE);
            }
          }
        }
      }

      // Draw particles
      if (settingsRef.current.showParticles) {
        particles.forEach((particle) => {
          const speed = Math.sqrt(
            particle.vx * particle.vx + particle.vy * particle.vy
          );
          const hue = 60 + speed * 30;
          const alpha = particle.life * 0.8;

          ctx.fillStyle = `hsla(${hue}, 80%, 60%, ${alpha})`;
          ctx.beginPath();
          ctx.arc(particle.x, particle.y, 1.5, 0, Math.PI * 2);
          ctx.fill();
        });
      }
    }

    // Draw obstacles (same for both modes)
    ctx.fillStyle =
      drawingMode === "draw"
        ? "#ff6b6b"
        : visualizationMode === "pressure"
        ? "#222"
        : "#444";
    for (let i = 0; i < GRID_WIDTH; i++) {
      for (let j = 0; j < GRID_HEIGHT; j++) {
        const idx = j * GRID_WIDTH + i;
        if (obstacles[idx]) {
          ctx.fillRect(i * CELL_SIZE, j * CELL_SIZE, CELL_SIZE, CELL_SIZE);
        }
      }
    }

    // Draw bezier curve preview
    if (drawingMode === "bezier" && currentBezier.length > 0) {
      ctx.strokeStyle = "#00ff00";
      ctx.lineWidth = 2;
      ctx.setLineDash([5, 5]);

      // Draw control points
      currentBezier.forEach((point, index) => {
        ctx.fillStyle = index === 0 || index === 3 ? "#00ff00" : "#ffff00";
        ctx.beginPath();
        ctx.arc(point.x, point.y, 4, 0, Math.PI * 2);
        ctx.fill();
      });

      // Draw preview curve if we have enough points
      if (currentBezier.length >= 2) {
        ctx.beginPath();
        ctx.moveTo(currentBezier[0].x, currentBezier[0].y);
        for (let i = 1; i < currentBezier.length; i++) {
          ctx.lineTo(currentBezier[i].x, currentBezier[i].y);
        }
        ctx.stroke();
      }

      ctx.setLineDash([]);
    }
  }, []); // Remove dependencies - render function will use current state via closure

  // Update particles
  const updateParticles = () => {
    const particles = simulationState.current.particles;
    const dt = 0.016;

    // Adjust particle count dynamically
    const targetCount = settingsRef.current.particleCount;
    const currentCount = particles.length;

    if (currentCount < targetCount) {
      // Add particles
      for (let i = currentCount; i < targetCount; i++) {
        particles.push({
          x: Math.random() * GRID_WIDTH * CELL_SIZE,
          y: Math.random() * GRID_HEIGHT * CELL_SIZE,
          vx: 0,
          vy: 0,
          life: 1.0,
        });
      }
    } else if (currentCount > targetCount) {
      // Remove particles
      particles.splice(targetCount);
    }

    particles.forEach((particle) => {
      // Get velocity at particle position
      const gridX = particle.x / CELL_SIZE;
      const gridY = particle.y / CELL_SIZE;

      if (
        gridX >= 0 &&
        gridX < GRID_WIDTH &&
        gridY >= 0 &&
        gridY < GRID_HEIGHT
      ) {
        const vx = interpolateVelocity(
          gridX,
          gridY,
          simulationState.current.velocityX
        );
        const vy = interpolateVelocity(
          gridX,
          gridY,
          simulationState.current.velocityY
        );

        particle.vx = vx;
        particle.vy = vy;
        particle.x += vx * dt * CELL_SIZE;
        particle.y += vy * dt * CELL_SIZE;
      }

      // Wrap around or respawn particles
      if (
        particle.x > GRID_WIDTH * CELL_SIZE ||
        particle.x < 0 ||
        particle.y > GRID_HEIGHT * CELL_SIZE ||
        particle.y < 0
      ) {
        particle.x = Math.random() * CELL_SIZE * 5;
        particle.y = Math.random() * GRID_HEIGHT * CELL_SIZE;
        particle.life = 1.0;
      }

      particle.life *= 0.999;
      if (particle.life < 0.1) {
        particle.life = 1.0;
      }
    });
  };

  // Drawing functions
  const getCanvasCoordinates = (e) => {
    const canvas = canvasRef.current;
    if (!canvas) return { x: 0, y: 0 };

    const rect = canvas.getBoundingClientRect();
    const scaleX = canvas.width / rect.width;
    const scaleY = canvas.height / rect.height;

    return {
      x: (e.clientX - rect.left) * scaleX,
      y: (e.clientY - rect.top) * scaleY,
    };
  };

  const drawObstacle = (x, y, isErasing = false) => {
    const obstacles = simulationState.current.obstacles;
    const gridX = Math.floor(x / CELL_SIZE);
    const gridY = Math.floor(y / CELL_SIZE);
    const radius = Math.max(1, Math.floor(brushSize / (CELL_SIZE * 2)));

    // Ensure we're within bounds
    if (gridX < 0 || gridX >= GRID_WIDTH || gridY < 0 || gridY >= GRID_HEIGHT) {
      return;
    }

    for (
      let i = Math.max(0, gridX - radius);
      i <= Math.min(GRID_WIDTH - 1, gridX + radius);
      i++
    ) {
      for (
        let j = Math.max(0, gridY - radius);
        j <= Math.min(GRID_HEIGHT - 1, gridY + radius);
        j++
      ) {
        const dx = i - gridX;
        const dy = j - gridY;
        if (dx * dx + dy * dy <= radius * radius) {
          const idx = j * GRID_WIDTH + i;
          obstacles[idx] = !isErasing;

          // Clear velocity at obstacle locations
          if (!isErasing) {
            simulationState.current.velocityX[idx] = 0;
            simulationState.current.velocityY[idx] = 0;
          }
        }
      }
    }
  };

  // Bezier curve functions
  const drawBezierCurve = (points, thickness = 3) => {
    if (points.length < 4) return;

    const obstacles = simulationState.current.obstacles;

    // Draw bezier curve using De Casteljau's algorithm
    for (let t = 0; t <= 1; t += 0.01) {
      const point = getBezierPoint(points, t);
      const gridX = Math.floor(point.x / CELL_SIZE);
      const gridY = Math.floor(point.y / CELL_SIZE);

      // Draw thick line
      const radius = Math.floor(thickness / 2);
      for (
        let i = Math.max(0, gridX - radius);
        i <= Math.min(GRID_WIDTH - 1, gridX + radius);
        i++
      ) {
        for (
          let j = Math.max(0, gridY - radius);
          j <= Math.min(GRID_HEIGHT - 1, gridY + radius);
          j++
        ) {
          const dx = i - gridX;
          const dy = j - gridY;
          if (dx * dx + dy * dy <= radius * radius) {
            const idx = j * GRID_WIDTH + i;
            obstacles[idx] = true;
            simulationState.current.velocityX[idx] = 0;
            simulationState.current.velocityY[idx] = 0;
          }
        }
      }
    }
  };

  const getBezierPoint = (points, t) => {
    if (points.length === 4) {
      // Cubic bezier
      const [p0, p1, p2, p3] = points;
      const u = 1 - t;
      return {
        x:
          u * u * u * p0.x +
          3 * u * u * t * p1.x +
          3 * u * t * t * p2.x +
          t * t * t * p3.x,
        y:
          u * u * u * p0.y +
          3 * u * u * t * p1.y +
          3 * u * t * t * p2.y +
          t * t * t * p3.y,
      };
    }
    return points[0] || { x: 0, y: 0 };
  };

  const handleMouseDown = useCallback(
    (e) => {
      e.preventDefault();
      const coords = getCanvasCoordinates(e);

      if (drawingMode === "bezier") {
        const newPoint = { x: coords.x, y: coords.y };
        setCurrentBezier((prev) => {
          const updated = [...prev, newPoint];
          if (updated.length === 4) {
            // Complete the bezier curve
            drawBezierCurve(updated, brushSize / 2);
            setBezierPoints((prevBezier) => [...prevBezier, updated]);
            if (!isRunning) render();
            return [];
          }
          return updated;
        });
      } else if (drawingMode !== "none") {
        setIsDrawing(true);
        drawObstacle(coords.x, coords.y, drawingMode === "erase");
        if (!isRunning) {
          render();
        }
      }
    },
    [drawingMode, isRunning, brushSize, render]
  );

  const handleMouseMove = useCallback(
    (e) => {
      e.preventDefault();
      if (isDrawing && drawingMode !== "none" && drawingMode !== "bezier") {
        const coords = getCanvasCoordinates(e);
        drawObstacle(coords.x, coords.y, drawingMode === "erase");
        if (!isRunning) {
          render();
        }
      }
    },
    [isDrawing, drawingMode, isRunning, brushSize, render]
  );

  const handleMouseUp = useCallback((e) => {
    e.preventDefault();
    setIsDrawing(false);
  }, []);

  const handleMouseLeave = useCallback((e) => {
    e.preventDefault();
    setIsDrawing(false);
  }, []);

  const clearObstacles = () => {
    simulationState.current.obstacles.fill(false);
    setBezierPoints([]);
    setCurrentBezier([]);
    if (!isRunning) {
      render();
    }
  };
  const simulationStep = useCallback(() => {
    applyBoundaryConditions();
    applyViscosity();
    projectVelocity();
    advectVelocity();
    applyBoundaryConditions();
    updateParticles();
  }, []); // Remove settings dependency - functions will use current settings via closure

  // Animation loop
  const animate = useCallback(() => {
    if (isRunning) {
      simulationStep();
      render();
      animationRef.current = requestAnimationFrame(animate);
    }
  }, [isRunning, simulationStep, render]);

  // Initialize simulation once
  useEffect(() => {
    initializeObstacles();
    initializeParticles();
    render();
  }, [initializeObstacles, initializeParticles, render]); // Only run once on mount

  // Reinitialize only when quality changes (grid size changes)
  useEffect(() => {
    setIsRunning(false);
    setTimeout(() => {
      reinitializeSimulation();
      initializeObstacles();
      initializeParticles();
      render();
    }, 100);
  }, [
    initializeObstacles,
    initializeParticles,
    quality,
    reinitializeSimulation,
    render,
  ]); // Only depend on quality

  // Re-render when drawing mode, visualization, or settings change (without restarting simulation)
  useEffect(() => {
    if (!isRunning) {
      render();
    }
  }, [drawingMode, visualizationMode, settings, isRunning]);

  // Start/stop animation
  useEffect(() => {
    if (isRunning) {
      animate();
    } else {
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    }
    return () => {
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    };
  }, [isRunning, animate]);

  const handlePlayPause = () => {
    setIsRunning(!isRunning);
  };

  const handleReset = () => {
    setIsRunning(false);
    setTimeout(() => {
      initializeParticles();
      simulationState.current.velocityX.fill(0);
      simulationState.current.velocityY.fill(0);
      simulationState.current.pressure.fill(0);
      render();
    }, 100);
  };

  return (
    <div className="min-h-screen bg-gray-900 text-white p-4">
      <div className="max-w-7xl mx-auto">
        <div className="text-center mb-6">
          <h1 className="text-4xl font-bold mb-2 bg-gradient-to-r from-blue-400 to-cyan-400 bg-clip-text text-transparent">
            CFD Fluid Dynamics Simulator
          </h1>
          <p className="text-gray-400">
            Interactive computational fluid dynamics visualization
          </p>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
          {/* Main Canvas */}
          <div className="lg:col-span-3">
            <div className="bg-gray-800 rounded-lg p-4">
              <div className="flex items-center justify-between mb-4">
                <div>
                  <h2 className="text-xl font-semibold flex items-center gap-2">
                    <Wind className="w-5 h-5" />
                    Flow Visualization
                  </h2>
                  <p className="text-sm text-gray-400">
                    Grid: {GRID_WIDTH}×{GRID_HEIGHT} ({quality} quality)
                  </p>
                </div>
                <div className="flex gap-2 flex-wrap">
                  <div className="flex gap-1 bg-gray-700 rounded-lg p-1">
                    <button
                      onClick={() => setVisualizationMode("standard")}
                      className={`px-3 py-1 rounded text-sm transition-colors ${
                        visualizationMode === "standard"
                          ? "bg-blue-600 text-white"
                          : "text-gray-300 hover:text-white"
                      }`}
                    >
                      Standard
                    </button>
                    <button
                      onClick={() => setVisualizationMode("pressure")}
                      className={`px-3 py-1 rounded text-sm transition-colors ${
                        visualizationMode === "pressure"
                          ? "bg-purple-600 text-white"
                          : "text-gray-300 hover:text-white"
                      }`}
                    >
                      Pressure
                    </button>
                  </div>

                  <button
                    onClick={() =>
                      setDrawingMode(drawingMode === "draw" ? "none" : "draw")
                    }
                    className={`px-3 py-2 rounded-lg flex items-center gap-2 transition-colors ${
                      drawingMode === "draw"
                        ? "bg-red-600 hover:bg-red-700"
                        : "bg-gray-600 hover:bg-gray-700"
                    }`}
                  >
                    <Brush className="w-4 h-4" />
                    Draw
                  </button>
                  <button
                    onClick={() =>
                      setDrawingMode(
                        drawingMode === "bezier" ? "none" : "bezier"
                      )
                    }
                    className={`px-3 py-2 rounded-lg flex items-center gap-2 transition-colors ${
                      drawingMode === "bezier"
                        ? "bg-green-600 hover:bg-green-700"
                        : "bg-gray-600 hover:bg-gray-700"
                    }`}
                  >
                    <svg
                      className="w-4 h-4"
                      viewBox="0 0 24 24"
                      fill="none"
                      stroke="currentColor"
                    >
                      <path d="M3 12c0 0 3-6 9-6s9 6 9 6-3 6-9 6-9-6-9-6z" />
                    </svg>
                    Bezier
                  </button>
                  <button
                    onClick={() =>
                      setDrawingMode(drawingMode === "erase" ? "none" : "erase")
                    }
                    className={`px-3 py-2 rounded-lg flex items-center gap-2 transition-colors ${
                      drawingMode === "erase"
                        ? "bg-cyan-600 hover:bg-cyan-700"
                        : "bg-gray-600 hover:bg-gray-700"
                    }`}
                  >
                    <Eraser className="w-4 h-4" />
                    Erase
                  </button>
                  {drawingMode === "bezier" && currentBezier.length > 0 && (
                    <button
                      onClick={() => setCurrentBezier([])}
                      className="px-3 py-2 bg-yellow-600 hover:bg-yellow-700 rounded-lg transition-colors text-sm"
                    >
                      Cancel Curve
                    </button>
                  )}
                  <button
                    onClick={clearObstacles}
                    className="px-3 py-2 bg-orange-600 hover:bg-orange-700 rounded-lg transition-colors text-sm"
                  >
                    Clear All
                  </button>
                  <button
                    onClick={handlePlayPause}
                    className={`px-4 py-2 rounded-lg flex items-center gap-2 ${
                      isRunning
                        ? "bg-red-600 hover:bg-red-700"
                        : "bg-green-600 hover:bg-green-700"
                    } transition-colors`}
                  >
                    {isRunning ? (
                      <Pause className="w-4 h-4" />
                    ) : (
                      <Play className="w-4 h-4" />
                    )}
                    {isRunning ? "Pause" : "Play"}
                  </button>
                  <button
                    onClick={handleReset}
                    className="px-4 py-2 bg-gray-600 hover:bg-gray-700 rounded-lg flex items-center gap-2 transition-colors"
                  >
                    <RotateCcw className="w-4 h-4" />
                    Reset
                  </button>
                </div>
              </div>
              <div className="bg-black rounded-lg overflow-hidden relative">
                <canvas
                  ref={canvasRef}
                  width={GRID_WIDTH * CELL_SIZE}
                  height={GRID_HEIGHT * CELL_SIZE}
                  className={`w-full h-auto ${
                    quality === "ultra"
                      ? "max-h-[700px]"
                      : quality === "high"
                      ? "max-h-[650px]"
                      : "max-h-[600px]"
                  } ${
                    drawingMode === "draw"
                      ? "cursor-crosshair"
                      : drawingMode === "erase"
                      ? "cursor-pointer"
                      : "cursor-default"
                  }`}
                  onMouseDown={handleMouseDown}
                  onMouseMove={handleMouseMove}
                  onMouseUp={handleMouseUp}
                  onMouseLeave={handleMouseLeave}
                />
                {drawingMode !== "none" && (
                  <div className="absolute top-2 left-2 bg-black bg-opacity-75 text-white px-2 py-1 rounded text-sm">
                    {drawingMode === "draw" && "🖌️ Drawing obstacles"}
                    {drawingMode === "erase" && "🧽 Erasing obstacles"}
                    {drawingMode === "bezier" && (
                      <div>
                        📐 Bezier curve: Click 4 points
                        {currentBezier.length > 0 && (
                          <div className="text-xs">
                            Point {currentBezier.length}/4 -{" "}
                            {currentBezier.length === 1
                              ? "Start point"
                              : currentBezier.length === 2
                              ? "Control point 1"
                              : currentBezier.length === 3
                              ? "Control point 2"
                              : "End point"}
                          </div>
                        )}
                      </div>
                    )}
                  </div>
                )}
              </div>
            </div>
          </div>

          {/* Controls */}
          <div className="space-y-4">
            <div className="bg-gray-800 rounded-lg p-4">
              <h3 className="text-lg font-semibold mb-4 flex items-center gap-2">
                <Settings className="w-5 h-5" />
                Simulation Parameters
              </h3>

              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium mb-2">
                    Simulation Quality
                  </label>
                  <select
                    value={quality}
                    onChange={(e) => setQuality(e.target.value)}
                    className="w-full bg-gray-700 border border-gray-600 rounded px-3 py-2 text-white"
                  >
                    <option value="low">Low (80×40) - Fast</option>
                    <option value="medium">Medium (120×60) - Balanced</option>
                    <option value="high">High (160×80) - Detailed</option>
                    <option value="ultra">Ultra (200×100) - Maximum</option>
                  </select>
                  <p className="text-xs text-gray-500 mt-1">
                    Higher quality = more accurate simulation but slower
                    performance
                  </p>
                </div>

                <div>
                  <label className="block text-sm font-medium mb-2">
                    {drawingMode === "bezier"
                      ? "Curve Thickness"
                      : "Brush Size"}
                    : {brushSize}px
                  </label>
                  <input
                    type="range"
                    min="5"
                    max="50"
                    step="5"
                    value={brushSize}
                    onChange={(e) => setBrushSize(parseInt(e.target.value))}
                    className="w-full"
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium mb-2">
                    Arrow Size: {settings.arrowSize.toFixed(1)}x
                  </label>
                  <input
                    type="range"
                    min="0.5"
                    max="3.0"
                    step="0.1"
                    value={settings.arrowSize}
                    onChange={(e) =>
                      setSettings((prev) => ({
                        ...prev,
                        arrowSize: parseFloat(e.target.value),
                      }))
                    }
                    className="w-full"
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium mb-2">
                    Wind Speed: {settings.windSpeed.toFixed(1)}
                  </label>
                  <input
                    type="range"
                    min="0"
                    max="15"
                    step="0.5"
                    value={settings.windSpeed}
                    onChange={(e) =>
                      setSettings((prev) => ({
                        ...prev,
                        windSpeed: parseFloat(e.target.value),
                      }))
                    }
                    className="w-full"
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium mb-2">
                    Viscosity: {settings.viscosity.toFixed(3)}
                  </label>
                  <input
                    type="range"
                    min="0"
                    max="0.1"
                    step="0.005"
                    value={settings.viscosity}
                    onChange={(e) =>
                      setSettings((prev) => ({
                        ...prev,
                        viscosity: parseFloat(e.target.value),
                      }))
                    }
                    className="w-full"
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium mb-2">
                    Particle Count: {settings.particleCount}
                  </label>
                  <input
                    type="range"
                    min="100"
                    max="20000"
                    step="100"
                    value={settings.particleCount}
                    onChange={(e) =>
                      setSettings((prev) => ({
                        ...prev,
                        particleCount: parseInt(e.target.value),
                      }))
                    }
                    className="w-full"
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium mb-2">
                    Default Obstacle Size: {settings.obstacleSize}
                  </label>
                  <input
                    type="range"
                    min="0"
                    max="100"
                    step="10"
                    value={settings.obstacleSize}
                    onChange={(e) =>
                      setSettings((prev) => ({
                        ...prev,
                        obstacleSize: parseInt(e.target.value),
                      }))
                    }
                    className="w-full"
                  />
                  <p className="text-xs text-gray-500 mt-1">
                    Set to 0 to start with no default obstacle
                  </p>
                </div>
              </div>
            </div>

            <div className="bg-gray-800 rounded-lg p-4">
              <h3 className="text-lg font-semibold mb-4 flex items-center gap-2">
                <Droplets className="w-5 h-5" />
                Visualization
              </h3>

              <div className="space-y-3">
                <div className="mb-3">
                  <div className="text-sm font-medium mb-2">
                    Visualization Mode
                  </div>
                  <div className="text-xs text-gray-400 mb-2">
                    {visualizationMode === "standard"
                      ? "Standard flow visualization with particles and velocity vectors"
                      : "Enhanced pressure field visualization with contour lines"}
                  </div>
                </div>

                {visualizationMode === "standard" && (
                  <>
                    <label className="flex items-center gap-2">
                      <input
                        type="checkbox"
                        checked={settings.showVelocityField}
                        onChange={(e) =>
                          setSettings((prev) => ({
                            ...prev,
                            showVelocityField: e.target.checked,
                          }))
                        }
                        className="rounded"
                      />
                      <span className="text-sm">Velocity Field</span>
                    </label>

                    <label className="flex items-center gap-2">
                      <input
                        type="checkbox"
                        checked={settings.showPressureField}
                        onChange={(e) =>
                          setSettings((prev) => ({
                            ...prev,
                            showPressureField: e.target.checked,
                          }))
                        }
                        className="rounded"
                      />
                      <span className="text-sm">Pressure Field</span>
                    </label>

                    <label className="flex items-center gap-2">
                      <input
                        type="checkbox"
                        checked={settings.showParticles}
                        onChange={(e) =>
                          setSettings((prev) => ({
                            ...prev,
                            showParticles: e.target.checked,
                          }))
                        }
                        className="rounded"
                      />
                      <span className="text-sm">Particle Traces</span>
                    </label>
                  </>
                )}

                {visualizationMode === "pressure" && (
                  <label className="flex items-center gap-2">
                    <input
                      type="checkbox"
                      checked={settings.showVelocityField}
                      onChange={(e) =>
                        setSettings((prev) => ({
                          ...prev,
                          showVelocityField: e.target.checked,
                        }))
                      }
                      className="rounded"
                    />
                    <span className="text-sm">Velocity Vectors</span>
                  </label>
                )}
              </div>
            </div>

            <div className="bg-gray-800 rounded-lg p-4">
              <h3 className="text-lg font-semibold mb-3">How to Use</h3>
              <div className="text-sm text-gray-400 space-y-2">
                <p>
                  <strong>� VDisualization:</strong> Switch between Standard and
                  Pressure modes
                </p>
                <p>
                  <strong>🖌️ Draw Mode:</strong> Click and drag to draw
                  obstacles
                </p>
                <p>
                  <strong>🧽 Erase Mode:</strong> Click and drag to remove
                  obstacles
                </p>
                <p>
                  <strong>▶️ Simulation:</strong> Press Play to start fluid
                  simulation
                </p>
                <p>
                  <strong>⚙️ Parameters:</strong> Adjust settings to see
                  different flow behaviors
                </p>
              </div>
            </div>

            <div className="bg-gray-800 rounded-lg p-4">
              <h3 className="text-lg font-semibold mb-3">Pressure Mode</h3>
              <div className="text-sm text-gray-400 space-y-2">
                <div className="flex items-center gap-2">
                  <div className="w-3 h-3 bg-red-500 rounded"></div>
                  <span>High Pressure (Red/Yellow)</span>
                </div>
                <div className="flex items-center gap-2">
                  <div className="w-3 h-3 bg-blue-500 rounded"></div>
                  <span>Low Pressure (Blue/Cyan)</span>
                </div>
                <p className="text-xs">
                  White dots show pressure contour lines
                </p>
              </div>
            </div>

            {quality === "ultra" && (
              <div className="bg-yellow-900 border border-yellow-600 rounded-lg p-4">
                <h3 className="text-lg font-semibold mb-2 text-yellow-300">
                  ⚠️ Performance Notice
                </h3>
                <p className="text-sm text-yellow-200">
                  Ultra quality uses a 200×100 grid which may impact performance
                  on slower devices. Consider using High quality for the best
                  balance of detail and performance.
                </p>
              </div>
            )}

            <div className="bg-gray-800 rounded-lg p-4">
              <h3 className="text-lg font-semibold mb-3">About</h3>
              <p className="text-sm text-gray-400">
                This CFD simulator demonstrates fluid flow around custom
                obstacles using the Navier-Stokes equations. Draw your own
                shapes to see how different geometries affect flow patterns,
                vortex formation, and pressure distribution.
              </p>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default CFDSimulator;
