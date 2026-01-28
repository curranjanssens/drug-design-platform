"""
Live Narrator Service

Uses a fast LLM (Haiku) to generate real-time, human-readable summaries
of what's happening during drug design. Instead of generic status messages,
this provides actual insight into the work being done.

The narrator receives chunks of actual work (prompts, responses, data)
and distills them into concise, informative updates for the user.
"""

import asyncio
import logging
import os
import time
from typing import Optional, Callable, List, Dict, Any
from dataclasses import dataclass, field
from collections import deque

import httpx

logger = logging.getLogger(__name__)


@dataclass
class WorkChunk:
    """A piece of work to be narrated."""
    timestamp: float
    phase: str
    content: str  # The actual work content (prompt, response, data)
    context: Dict[str, Any] = field(default_factory=dict)


class LiveNarrator:
    """
    Generates live narrative updates about ongoing work.

    Uses a fast model (Haiku) to summarize actual work into
    user-friendly updates every few seconds.
    """

    def __init__(self, api_key: str = None):
        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        self._work_buffer: deque = deque(maxlen=10)  # Recent work chunks
        self._last_narration_time = 0
        self._min_interval = 8  # Minimum seconds between narrations
        self._callback: Optional[Callable] = None
        self._running = False
        self._task: Optional[asyncio.Task] = None

        # Context about the overall task
        self._target_name = ""
        self._mechanism = ""
        self._phase = ""

    def set_callback(self, callback: Callable):
        """Set the callback for narrative updates."""
        self._callback = callback

    def set_context(self, target: str = "", mechanism: str = "", phase: str = ""):
        """Update the context for narration."""
        if target:
            self._target_name = target
        if mechanism:
            self._mechanism = mechanism
        if phase:
            self._phase = phase

    def add_work(self, phase: str, content: str, context: Dict = None):
        """Add a chunk of work to be potentially narrated."""
        self._work_buffer.append(WorkChunk(
            timestamp=time.time(),
            phase=phase,
            content=content[:2000],  # Limit size
            context=context or {}
        ))
        self._phase = phase

    async def start(self):
        """Start the narrator background task."""
        self._running = True
        self._task = asyncio.create_task(self._narration_loop())

    async def stop(self):
        """Stop the narrator."""
        self._running = False
        if self._task:
            self._task.cancel()
            try:
                await self._task
            except asyncio.CancelledError:
                pass

    async def _narration_loop(self):
        """Background loop that generates narrations periodically."""
        while self._running:
            try:
                await asyncio.sleep(2)  # Check every 2 seconds

                # Check if enough time has passed and we have new work
                now = time.time()
                if now - self._last_narration_time < self._min_interval:
                    continue

                if not self._work_buffer:
                    continue

                # Get recent work that hasn't been narrated
                recent_work = list(self._work_buffer)
                if not recent_work:
                    continue

                # Generate narration
                narration = await self._generate_narration(recent_work)

                if narration and self._callback:
                    await self._callback(self._phase, narration, {
                        "type": "narrative",
                        "target": self._target_name
                    })

                self._last_narration_time = now
                self._work_buffer.clear()

            except asyncio.CancelledError:
                break
            except Exception as e:
                logger.warning(f"Narration error: {e}")

    async def _generate_narration(self, work_chunks: List[WorkChunk]) -> Optional[str]:
        """Generate a narrative summary of recent work."""
        if not self.api_key:
            return None

        # Build context from work chunks
        work_summary = []
        for chunk in work_chunks[-5:]:  # Last 5 chunks
            work_summary.append(f"[{chunk.phase}] {chunk.content[:500]}")

        work_text = "\n---\n".join(work_summary)

        prompt = f"""You are narrating a live drug design process to a scientist watching the UI.

TARGET: {self._target_name or "Being identified"}
MECHANISM: {self._mechanism or "Being analyzed"}
CURRENT PHASE: {self._phase}

RECENT WORK BEING DONE:
{work_text}

Write a single, concise sentence (max 100 chars) that tells the user what's ACTUALLY happening right now.
- Be specific about what was found/done (mention actual compounds, targets, scores if available)
- Use present tense
- Don't be generic - extract real information from the work
- If analyzing compounds, mention what kind
- If found bioactivity data, mention the values
- If validating molecules, mention how many passed/failed

Examples of GOOD narrations:
- "Found PF-04457845 (Ki=0.7nM) as lead reference compound"
- "Identified urea warhead targeting Ser241 - analyzing leaving group chemistry"
- "Generated 10 candidates, 7 passed topology validation"
- "Best candidate scores 83% - piperidine-biaryl scaffold with CF3"
- "ChEMBL returned 47 active compounds, analyzing SAR patterns"

Examples of BAD narrations (too generic):
- "Analyzing target..."
- "Processing compounds..."
- "Running design..."

Return ONLY the narration sentence, nothing else."""

        try:
            async with httpx.AsyncClient(timeout=10.0) as client:
                response = await client.post(
                    "https://api.anthropic.com/v1/messages",
                    headers={
                        "x-api-key": self.api_key,
                        "anthropic-version": "2023-06-01",
                        "content-type": "application/json"
                    },
                    json={
                        "model": "claude-3-5-haiku-20241022",  # Fast model for quick narrations
                        "max_tokens": 150,
                        "messages": [{"role": "user", "content": prompt}]
                    }
                )

                if response.status_code == 200:
                    data = response.json()
                    narration = data["content"][0]["text"].strip()
                    # Clean up - remove quotes if present
                    narration = narration.strip('"\'')
                    return narration[:150]  # Limit length
                else:
                    logger.warning(f"Narration API error: {response.status_code}")
                    return None

        except Exception as e:
            logger.warning(f"Narration generation failed: {e}")
            return None

    async def narrate_now(self, content: str, phase: str = None) -> Optional[str]:
        """Generate an immediate narration for specific content."""
        if phase:
            self._phase = phase

        # Create a single work chunk and narrate it
        chunk = WorkChunk(
            timestamp=time.time(),
            phase=phase or self._phase,
            content=content
        )

        narration = await self._generate_narration([chunk])

        if narration and self._callback:
            await self._callback(self._phase, narration, {
                "type": "narrative",
                "target": self._target_name
            })

        return narration


# Singleton instance
live_narrator = LiveNarrator()
