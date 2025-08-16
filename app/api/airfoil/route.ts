import { NextRequest, NextResponse } from "next/server";

export async function GET(request: NextRequest) {
  const { searchParams } = new URL(request.url);
  const airfoil = searchParams.get("airfoil");

  if (!airfoil) {
    return NextResponse.json(
      { error: "Airfoil parameter is required" },
      { status: 400 }
    );
  }

  try {
    const response = await fetch(
      `http://airfoiltools.com/airfoil/seligdatfile?airfoil=${airfoil}`,
      {
        headers: {
          "User-Agent": "Mozilla/5.0 (compatible; CFD-Simulator/1.0)",
        },
      }
    );

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    const data = await response.text();

    return new NextResponse(data, {
      headers: {
        "Content-Type": "text/plain",
        "Access-Control-Allow-Origin": "*",
        "Access-Control-Allow-Methods": "GET",
        "Access-Control-Allow-Headers": "Content-Type",
      },
    });
  } catch (error) {
    console.error("Error fetching airfoil data:", error);
    return NextResponse.json(
      { error: "Failed to fetch airfoil data" },
      { status: 500 }
    );
  }
}
